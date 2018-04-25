#!/usr/bin/env python2.7

import sys 
import argparse as ap

# help
if len(sys.argv) < 2:
	sys.argv.append("-h")
if sys.argv[1] in ["-h", "-help", "--help", "usage", "getopt"]:
	sys.exit('''

---
add-read-groups.py
version 1.0	April 24, 2018

Matteo Schiavinato
BOKU - Vienna University of Natural Resources and Life Sciences
---

Warning: this program will modify the SAM header and the read TAGS (fields from 12 on).
The modifications will be the same for the whole file (i.e. each read will get the same read group) 
It was originally intended to be used in GATK data preparation. 


Usage:
samtools view -h | add-read-groups.py [OPTIONS]

Mandatory:
-id	--rgid		ID				[-]
-lb	--rglb		Library				[lib1]
-pl	--rgpl		Platform 			[illumina]
-pu	--rgpu		Platform Unit			[unit1]
-sm	--rgsm		Sample				[-]

Optional:
-cn	--rgcn		Sequencing center		[-]
-ds	--rgds		Description			[-]
-dt	--rgdt		Run date			[-]
-pi	--rgpi		Predicted insert size		[-]
-pg	--rgpg		Program group			[-]
-pm	--rgpm		Platform model			[-]

''')

p = ap.ArgumentParser()
p.add_argument("-id", "--rgid")
p.add_argument("-lb", "--rglb", default="lib1")
p.add_argument("-pl", "--rgpl", default="illumina")
p.add_argument("-pu", "--rgpu", default="unit1")
p.add_argument("-sm", "--rgsm")
p.add_argument("-cn", "--rgcn")
p.add_argument("-ds", "--rgds")
p.add_argument("-dt", "--rgdt")
p.add_argument("-pi", "--rgpi")
p.add_argument("-pg", "--rgpg")
p.add_argument("-pm", "--rgpm")
args = p.parse_args()

Conversions = {
"rgid":"ID",
"rglb":"LB",
"rgpl":"PL",
"rgpu":"PU",
"rgsm":"SM",
"rgcn":"CN",
"rgds":"DS",
"rgdt":"DT",
"rgpi":"PI",
"rgpg":"PG",
"rgpm":"PM"}

new_RG = ["@RG"]
for arg in vars(args):
	value = getattr(args, arg)
	key = Conversions[str(arg)]
	if key == "ID":
		RGID=value
	if value != None:
		pair = str(key) + ":" + str(value)
		new_RG.append(pair)

new_RG_line = "\t".join(new_RG) + "\n"

counter = 0
status="not printed"
INPUT = sys.stdin
for line in sys.stdin:
	if line[0:3] == "@RG":
		status="printed"
		sys.stdout.write(new_RG_line)
	else:
		tags = line.rstrip("\b\r\n").split("\t")[11:]
		dic = {}
		for tag in tags:
			tag_key = ":".join([tag.split(":")[i] for i in (0,1)])
			tag_value = tag.split(":")[2]
			dic[tag_key] = tag_value
		if ("RG:Z" not in dic.keys()) or (dic["RG:Z"] != RGID):
			counter += 1
			dic["RG:Z"] = RGID
		newtags = "\t".join([str(tag_key)+":"+str(dic[tag_key]) for tag_key in dic])
		newline = "\t".join(line.rstrip("\b\r\n").split("\t")[:11]) + "\t" + newtags + "\n"
		sys.stdout.write(newline)

if status=="not_printed":
	sys.stderr.write('''

ERROR: Read group information was not added because an @RG line was not found.
See https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups

''')

sys.stderr.write('''

Read group updated in the header: DONE
Header read group line:
{0}

Read group modified for {1} mapping scores: DONE
Read group assigned:
{2}

'''.format(new_RG_line, counter, RGID))

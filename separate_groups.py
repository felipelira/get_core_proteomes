#!/usr/bin/env python

import subprocess
import sys

"""
List was created using this table:

GCF_000026185.1.gbff	P
GCF_000026985.1.gbff	P
GCF_000091565.1.gbff	P
GCF_000165815.1.gbff	C
GCF_000240705.2.gbff	C
GCF_000336255.1.gbff	C
GCF_000367545.1.gbff	C
GCF_000367585.1.gbff	C

if you can use it as an example
L = [["GCF_000026185.1.gbff","P"], ["GCF_000026985.1.gbff",	"P"], ["GCF_000091565.1.gbff","C"], ["GCF_000165815.1.gbff","P"], ["GCF_000240705.2.gbff","C"], ["GCF_000336255.1.gbff","C"], ["GCF_000367545.1.gbff","P"],["GCF_000367585.1.gbff","P"]]

"""

groups = {}
lists = []
input_parse = []
L = sys.argv[1]
	

for line in open(L, "r"):
	group = line.strip().split("\t")[1]
	genome = line.strip().split("\t")[0]
	if group not in groups:
		groups[group] = [genome]
		lists.append(group)
	else:
		groups[group].append(genome)


for g in lists:
	G = "list_" + g + ".txt"
	input_parse.append(G)
	list = open(G, "w")

	if g in groups.iterkeys():
		if g == "P":
			list.write('\n'.join(groups[g]))
		elif g == "C":
			list.write('\n'.join(groups[g]))

subprocess.call(["parse_pangenome_matrix.pl", "-m", "some_gbf_homologues_compared/pangenome_matrix_t0.tab", "-A", "list_P.txt", "-B", "list_C.txt", "-g"])




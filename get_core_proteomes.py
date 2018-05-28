#!/usr/bin/env python

"""
extract_plasmids_from_reports.py fragments a fasta sequence in small parts with fixed length
to be used in a blast search or bowtie alignment against a reference genome.

directory = directory used to calculate the core and pangenome clusters

		usage: get_core_proteomes.py [directory]
"""

import sys, os, errno
import subprocess
import argparse
from Bio import SeqIO
from optparse import OptionParser

#############  GLOBALS  #############

folders = []	# store the names of the directories previousl y created by 'get_homologues'
core_clusters = []	# store the names of files that correspond to the core genomes in core_clusters = []
species = set()

#############  FUNCTIONS  #############
# detect folders to be used by compare_clusters from "get_homologues.pl" package
def pipeline(inhandle):

	for d in os.listdir(os.path.join(inhandle)):
		element = os.path.join(inhandle, d)
		if os.path.isdir(element):
			if element.endswith("tmp"):
				pass
			else:
				folders.append(element)

	infolder = ",".join(folders)
	output = inhandle + "_compared"

	# run subprocesses to compare clusters and generate the list of genes that compose the core genome
	subprocess.call(["compare_clusters.pl", "-d"  , infolder, "-o", output, "-m", "-T"])
	matrix = os.path.join(output, "pangenome_matrix_t0.tab")			
	subprocess.call(["parse_pangenome_matrix.pl", "-m", matrix, "-s"])

	### EXTRACT CORE PROTEOMES ###
	# read directly the list of clusters of the core genome and store the name in list
	list_clusters = os.path.join(output, "pangenome_matrix_t0__core_list.txt")
	core = open(os.path.join(list_clusters), "r")
	for line in core:
		faa = line.strip()	# file with core genes
		cluster = os.path.join(output,faa)
		core_clusters.append(cluster)

	for faa in core_clusters:
		sequences = []
		if os.path.isfile(faa):
			parse = SeqIO.parse(faa, "fasta")
			for record in parse:
				sequences.append(record.id)
				desc = record.description
				split_desc = desc.strip().split("|")
				sp = split_desc[1].replace("[", "").replace("]", "")
				species.add(sp)
	# create folder to store proteomes
	proteomes = output + "_core_proteomes"
	try:
		os.mkdir(proteomes)
	except:
		pass

	for sp in species:
		proteome_core = open(os.path.join(proteomes,(sp.split("/")[0] + ".faa").replace(" ", "_")), "w")
		#out_fasta = open(os.path.join(proteome_core), "w")
		for faa in core_clusters:
			f = open(os.path.join(faa), "r")	# create file to store sequences
			parse = SeqIO.parse(faa, "fasta")

			for record in parse:
				d = record.description
				split_d = d.strip().split("|")
				s = split_d[1].replace("[", "").replace("]", "")
				#print faa, record.id
				if sp in s:
					seq = ">" + d + "\n" + str(record.seq) + "\n"
					proteome_core.write(seq)
	run_CVTree(proteomes)
	differential_genes()

# Run "CV" and "TREE" scripts
def run_CVTree(str_n):
	# name of the folder to store the list of proteomes and the outputs generated
	out_cvtree = str_n + "_CVTree"

	try:
		os.makedirs(out_cvtree)
	except OSError as e:
		if e.errno != errno.EEXIST:
			raise

	inList = os.path.join(out_cvtree,"list_CVTree.txt")

	with open(inList, "w") as inlist:
		for genome in os.listdir(str_n):
			inlist.write("%s\n" % (genome))
	inlist.close()

	print "File list_CVTree.txt created: ", inList

	# Calculate the matrix and dendrogram using core proteomes
	# Use two different lengths	
	kmers = ["7", "14"]
	for i in kmers:
		print "Running CVTree using k-mers length:", i, "\n"
		out_dist = os.path.join(out_cvtree ,("dist_kmer_" + i + ".dist"))
		out_tree = os.path.join(out_cvtree ,("tree_kmer_" + i + ".nwk"))
		tag = "cv" + i + ".gz"
		print "\tCalculating distances using kmer length: ", i, ""
		subprocess.call(["cv", "-I", str_n, "-i", inList, "-k",  i , "-O", out_cvtree]) # when files are proteomes
		print "\tGenerating tree using kmer length: ", i, ""
		subprocess.call(["tree", "-o", out_tree, "-I", out_cvtree, "-i", inList, "-s", tag, "-E", str_n, "-d", out_dist])

		
		
		
# detect genes present in the 
def differential_genes():

	groups = {}
	lists = []
	input_parse = []
	L = sys.argv[2]

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
		
		
		
# define the mode of pipeline, if it will run over a file or in a directory
def mode(inhandle):	
	if os.path.isdir(inhandle):
		pipeline(inhandle)
	else:
		argsCheck()

# check the number of arguments as input
def argsCheck():
	parser = argparse.ArgumentParser(prog='pipeline_pancore.py', usage='\n\t%(prog)s [directory]')
	parser.add_argument('-sys.argv[1]', help='Folder generated by get_homologues.pl.', metavar='')
	#parser.add_argument('-sys.argv[2]', help='List with group A', metavar='')
	#parser.add_argument('-sys.argv[3]', help='List with group B', metavar='')
	# TODO parser = argparse.ArgumentParser(prog='pipeline_pancore.py', usage='%(prog)s [directory] [list_of_groups]')
	parser.print_help()
	print "\n"

#############  CODE  #############
if __name__ == '__main__':
	if len(sys.argv) == 3: # script.py [directory] [list_of_groups]
		mode(sys.argv[1])
		
		#print len(sys.argv)

	elif len(sys.argv) != 3:
		print len(sys.argv)
		argsCheck()

# final project for MCB185 
# Compare codon usage by Kullback-Leibler distance
# K-L distance = sum pi*log2*(pi/qi)
# compare every gene to average codon usage in the genome

import argparse
import biotools
import re
import math

parser = argparse.ArgumentParser(
	description='Horiztonal gene transfer detector.')
parser.add_argument('--file', required=True, type=str,
	metavar='<str>', help='FASTA file')
parser.add_argument('--pseudo', required=False, type=float, default=1.0,
	metavar='<float>', help='pseudocount [%(default)f]')
arg = parser.parse_args()

def kl_distance(p, q):
	d = 0
	for i in p:
		d += p[i]*math.log2(p[i]/q[i])
	return d
	

#part 1 full codon table
total_counts = {}
sum = 0
for name, seq in biotools.read_fasta(arg.file):
	for i in range(0,len(seq), 3):
		codon = seq[i:i+3]
		sum += 1
		if codon in total_counts: total_counts[codon] += 1
		else:					  total_counts[codon]  = 1

freq_codons = {}	
for codon in total_counts:
	freq_codons[codon] = total_counts[codon]/sum


#part 2 compare every gene's table to the full codon table
for name, seq in biotools.read_fasta(arg.file):
	sum = 0
	gene_counts = {}
	gene_freq = {}
	for codon in total_counts:
		gene_counts[codon] = arg.pseudo
		sum += arg.pseudo
	for i in range(0, len(seq), 3):
		codon = seq[i:i+3]
		sum += 1
		gene_counts[codon] += 1
	for codon in gene_counts:
		gene_freq[codon] = gene_counts[codon]/sum	
	d = kl_distance(freq_codons, gene_freq)
	match = re.search('locus_tag=(\w+)', name)
	locus = match[1]
	print(locus, d)

# the higher the distance, the more likely it could've come from HGT


	

	
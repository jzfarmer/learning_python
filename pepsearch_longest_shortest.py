#!/usr/bin/env python3

import gzip
import sys
import math
# Write a program that finds peptidies within protein sequences
# Command line:
#	python3 pepsearch.py IAN
def read_fasta(filename):
	name = None
	seqs = []
	
	fp = None
	if filename == '-':
		fp = sys.stdin
	elif filename.endswith('.gz'):
		fp = gzip.open(filename, 'rt')
	else:
		fp = open(filename)

	for line in fp.readlines():
		line = line.rstrip()
		if line.startswith('>'):
			if len(seqs) > 0:
				seq = ''.join(seqs)
				yield(name, seq)
				name = line[1:]
				seqs = []
			else:
				name = line[1:]
		else:
			seqs.append(line)
	yield(name, ''.join(seqs))
	fp.close()

filename = sys.argv[1]
longest = -1
shortest = math.inf

for name, seq in read_fasta(sys.argv[1]):
	if len(seq) > longest : longest = len(seq)
	if len(seq) < shortest : shortest = len(seq)
print(shortest, longest)

"""
python3 pepsearch.py proteins.fasta.gz IAN | wc -w
	43
"""

#!/usr/bin/env python3

import sys
import gzip
import random

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

def gc(seq):
	count = 0
	for nt in seq:
		if nt == 'G' or nt == 'C':
			count += 1
	return count / len(seq)

def randseq(l, gc):
	dna = []
	for i in range(l):
		r = random.random()
		if r < gc:
			r = random.random()
			if r < 0.5: dna.append('G')
			else: dna.append('C')
		else:
			r = random.random()
			if r < 0.5: dna.append('T')
			else: dna.append('A')
	return ''.join(dna)

def kdh(aa):
	if   aa == 'A': return 1.8
	elif aa == 'C': return 2.5
	elif aa == 'D': return -3.5
	elif aa == 'E': return -3.5
	elif aa == 'F': return 2.8
	elif aa == 'G': return -0.4
	elif aa == 'H': return -3.2
	elif aa == 'I': return 4.5
	elif aa == 'K': return -3.9
	elif aa == 'L': return 3.8
	elif aa == 'M': return 1.9
	elif aa == 'N': return -3.5
	elif aa == 'P': return -1.6
	elif aa == 'Q': return -3.5
	elif aa == 'R': return -4.5
	elif aa == 'S': return -0.8
	elif aa == 'T': return -0.7
	elif aa == 'V': return 4.2
	elif aa == 'W': return -0.9
	elif aa == 'Y': return -1.3
	else: return 0

def skew(seq):
	g = 0
	c = 0
	for nt in seq:
		if nt == 'G':
			g += 1
		elif nt == 'C':
			c +=1
	return (g - c)/(g + c)
	
	
	
	
	
	
	
	
	
	
	
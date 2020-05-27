#!/usr/bin/env python3

import sys
import gzip
import random
import math

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
	
# percent gc content
def gc(seq):
	count = 0
	for nt in seq:
		if nt == 'G' or nt == 'C':
			count += 1
	return count / len(seq)

# generation of a random sequence
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

# kdh values
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
	
# gc skew of a sequence
def skew(seq):
	g = 0
	c = 0
	for nt in seq:
		if nt == 'G':
			g += 1
		elif nt == 'C':
			c +=1
	return (g - c)/(g + c)
	
# shannon entropy	
def entropy(seq):
	prob = [0]*4
	for nt in seq:
		if nt == 'A': prob[0] += 1
		elif nt == 'C': prob[1] += 1
		elif nt == 'G': prob[2] += 1
		elif nt == 'T': prob[3] += 1
		else: return None		
	h = 0
	for p in prob:
		if p == 0: h += 0
		else:
			h += (-p/len(seq)) * math.log2(p/len(seq))	
	return h	

# amino acid dictionary
gcode = {
	'AAA' : 'K',	'AAC' : 'N',	'AAG' : 'K',	'AAT' : 'N',
	'ACA' : 'T',	'ACC' : 'T',	'ACG' : 'T',	'ACT' : 'T',
	'AGA' : 'R',	'AGC' : 'S',	'AGG' : 'R',	'AGT' : 'S',
	'ATA' : 'I',	'ATC' : 'I',	'ATG' : 'M',	'ATT' : 'I',
	'CAA' : 'Q',	'CAC' : 'H',	'CAG' : 'Q',	'CAT' : 'H',
	'CCA' : 'P',	'CCC' : 'P',	'CCG' : 'P',	'CCT' : 'P',
	'CGA' : 'R',	'CGC' : 'R',	'CGG' : 'R',	'CGT' : 'R',
	'CTA' : 'L',	'CTC' : 'L',	'CTG' : 'L',	'CTT' : 'L',
	'GAA' : 'E',	'GAC' : 'D',	'GAG' : 'E',	'GAT' : 'D',
	'GCA' : 'A',	'GCC' : 'A',	'GCG' : 'A',	'GCT' : 'A',
	'GGA' : 'G',	'GGC' : 'G',	'GGG' : 'G',	'GGT' : 'G',
	'GTA' : 'V',	'GTC' : 'V',	'GTG' : 'V',	'GTT' : 'V',
	'TAA' : '*',	'TAC' : 'Y',	'TAG' : '*',	'TAT' : 'Y',
	'TCA' : 'S',	'TCC' : 'S',	'TCG' : 'S',	'TCT' : 'S',
	'TGA' : '*',	'TGC' : 'C',	'TGG' : 'W',	'TGT' : 'C',
	'TTA' : 'L',	'TTC' : 'F',	'TTG' : 'L',	'TTT' : 'F',
}

kdscale = { 'A': 1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C': 2.5,
       'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I': 4.5,
       'L': 3.8,'K':-3.9,'M': 1.9,'F': 2.8,'P':-1.6,
       'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V': 4.2 }

isscale = { 'A': 0.17,'R':0.81,'N':0.42,'D':1.23,'C':-0.24,
       'Q':0.58,'E':2.02,'G':2.02,'H':0.96,'I':-0.31,
       'L': -0.56,'K': 0.99,'M': -0.23,'F': -1.13,'P': 0.45,
       'S':0.13,'T':0.14,'W':-1.85,'Y':-0.94,'V': 0.07 }
       
osscale = { 'A': 0.50,'R':1.81,'N':0.85,'D':3.64,'C':-0.02,
       'Q':0.77,'E':3.63,'G':1.15,'H':2.33,'I':-1.12,
       'L': -1.25,'K': 2.80,'M': -0.67,'F': -1.71,'P': 0.14,
       'S':0.46,'T':0.25,'W':-2.09,'Y':-0.71,'V': 0.46 }

iosscale = { 'A': 0.33,'R':1.00,'N':0.43,'D':2.41,'C': 0.22,
       'Q':0.19,'E':1.61,'G':1.14,'H':1.37,'I':-0.81,
       'L': -0.69,'K': 1.81,'M': -0.44,'F': -0.58,'P': -0.31,
       'S':0.33,'T':0.11,'W':-0.24,'Y':0.23,'V': -0.53 }

ccscale = { 'A': -0.495,'R':4.383,'N':2.354,'D':9.573,'C':0.081,
       'Q':2.176,'E':3.173,'G':0.386,'H':2.029,'I':-0.528,
       'L': -0.342,'K': 2.101,'M': -0.324,'F': -0.370,'P': -0.322,
       'S':0.936,'T':0.853,'W':-0.270,'Y':1.677,'V': -0.308 }
       
       
def hydro(protein, method):
	if method == 'kd':
		scale = kdscale
	elif method == 'is':
		scale = isscale
	elif method == 'os':
		scale = osscale
	elif method == 'ios':
		scale = iosscale
	elif method == 'cc':
		scale = ccscale
	sum_score = 0
	for aa in protein:
		sum_score += scale[aa]
	return sum_score



	
	
	
	
	
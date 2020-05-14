#!/usr/bin/env python3

import argparse
import biotools as bt

# Write a program that predicts if a protein is trans-membrane
# Trans-membrane proteins have the following properties
#	Signal peptide: https://en.wikipedia.org/wiki/Signal_peptide
#	Hydrophobic regions(s): https://en.wikipedia.org/wiki/Transmembrane_protein
#	No prolines (alpha helix)
# Hydrophobicity is measued via Kyte-Dolittle
#	https://en.wikipedia.org/wiki/Hydrophilicity_plot
# For our purposes:
#	Signal peptide is 8 aa long, KD > 2.5, first 30 aa
#	Hydrophobic region is 11 aa long, KD > 2.0, after 30 aa
# setup
parser = argparse.ArgumentParser(
	description='Predicts trans-membrane proteins.')
# required arguments
parser.add_argument('--file', required=True, type=str,
	metavar='<path>', help='Protein file')
# optional arguments with default parameters
parser.add_argument('--win1', required=False, type=int, default=8,
	metavar='<int>', help='length of signal peptide [%(default)i]')
parser.add_argument('--win2', required=False, type=int, default=11,
	metavar='<int>', help='length of the transmembrane domain [%(default)i]')
parser.add_argument('--kd1', required=False, type=float, default=2.5,
	metavar='<float>', help='kd value for signal peptide  [%(default)f]')
parser.add_argument('--kd2', required=False, type=float, default=2.0,
	metavar='<float>', help='kd value for hydrophobicity [%(default)f]')
	
# finalization
arg = parser.parse_args()

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
	
def hydro(seq, w, t):
	for i in range(0, len(seq) -w + 1):
		pep = seq[i: i + w]
		sum = 0
		ave_hy = 0
		for aa in pep:
			sum += kdh(aa)
		ave_hy = sum / w
		if ave_hy > t and 'P' not in pep: return True
	return False	
		
for name, seq in bt.read_fasta(arg.file):
	nterm = seq[0:30]
	cterm = seq[30:-2]	
	if hydro(nterm, arg.kd1, arg.win1) and hydro(cterm, arg.kd2, arg.win2):
		print(name)

"""
python3 transmembrane.py proteins.fasta.gz
18w
Dtg
Krn
Lac
Mcr
PRY
Pxt
Pzl
QC
Ror
S1P
S2P
Spt
apn
bai
bdl
bou
bug
cue
drd
ft
grk
knk
ksh
m
nac
ort
rk
smo
thw
tsg
waw
zye
"""

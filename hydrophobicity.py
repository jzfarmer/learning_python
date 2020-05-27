#!/usr/bin/env python3
import argparse
import biotools as bt
# Write a program that computes hydrophobicity in a window
# Let the user choose the method (see below)
# https://en.wikipedia.org/wiki/Hydrophilicity_plot
# https://en.wikipedia.org/wiki/Hydrophobicity_scales

parser = argparse.ArgumentParser(
	description='calculate hydrophobicity')
# optional arguments with default parameters
parser.add_argument('--input', required=True, type=str,
	metavar='<path>', help='path to fasta file')
parser.add_argument('--window', required=False, type=int, default=15,
	metavar='<int>', help='integer argument for window size [%(default)i]')
parser.add_argument('--method', required=False, type=str, default='kd',
	metavar='<str>', help='Calculation method, kd, is, os, ios, cc [%(default)s]')
arg = parser.parse_args()

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
       
       
def hydro(protein, scale):
	sum_score = 0
	for aa in protein:
		sum_score += scale[aa]
	return sum_score

print(bt.hydro('ARNDC', 'kd'))
		

"""
python3 hydrophobicity.py --input proteins.fasta.gz --window 11 --method kd
"""

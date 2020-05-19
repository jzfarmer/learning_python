#!/usr/bin/env python3

import argparse
import biotools as bt
import sys

# Write a program that computes the amino acid composition of a protein file
# Use a dictionary
count = {}
total = 0

for id, protein in bt.read_fasta(sys.argv[1]):
	for aa in protein:
		if aa in count: count[aa] += 1
		else:           count[aa] = 1
	total += len(protein)
			
for aa in count:
	print(aa, count[aa]/total)
	


	
"""
python3 composition.py proteins.fasta.gz | sort -nk2
* 0.0017561333612916754
W 0.010255968606775905
C 0.019017913309169337
M 0.023765838900038944
H 0.027689991912051043
Y 0.02980558967138963
F 0.036174849474283316
N 0.04593281011293173
I 0.049422610310637154
D 0.052167270766557826
Q 0.05259413473923853
P 0.05463858850313034
T 0.05542491687385795
K 0.056080190516130966
G 0.05631234460653626
R 0.05732708264685618
V 0.05813962196327472
E 0.06519785519575833
A 0.07117020639247522
S 0.08295764311176347
L 0.09416843902585148
"""

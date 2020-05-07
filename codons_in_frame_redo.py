#!/usr/bin/env python3

# Print out all the codons for the sequence below in reading frame 1
# Use a 'for' loop

dna = 'ATAGCGAATATCTCTCATGAGAGGGAA'

for frame in range(3):
	print('frame')
	for i in range(0, len(dna)-2-frame, 3):
		print(i, dna[i + frame: i + 3 + frame])
	

"""
python3 codons.py
ATA
GCG
AAT
ATC
TCT
CAT
GAG
AGG
GAA
"""

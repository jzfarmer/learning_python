#!/usr/bin/env python3

# Print out all the codons for the sequence below in reading frame 1
# Use a 'for' loop

dna = 'ATAGCGAATATCTCTCATGAGAGGGAA'

for f in range(3):
	print ('frame',f)
	for i in range(f, len(dna) -2, 3):
		codon = dna[i:i+3]
		print(codon)
"""
frame0
ATA
GCG
AAT
ATC
TCT
CAT
GAG
AGG
GAA
frame1
TAG
CGA
ATA
....
frame2
"""

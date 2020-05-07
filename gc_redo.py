#!/usr/bin/env python3

# Write a program that computes the GC% of a DNA sequence
# Format the output for 2 decimal places
# Use all three formatting methods

dna = 'ACAGAGCCAGCAGATATACAGCAGATACTAT' # feel free to change

count = 0
for nt in dna:
	if nt == 'G' or nt == 'C': 
		count += 1
gc_content = count/len(dna)
print('%.2f' % (gc_content))



"""
python3 gc.py
0.42
0.42
0.42
"""

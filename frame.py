#!/usr/bin/env python3

# Write a program that prints out the position, frame, and letter of the DNA
# Try coding this with a single loop
# Try coding this with nested loops

dna = 'ATGGCCTTT'
for i in range (0,len(dna), 3):
	for j in range (0, 3):
		c = i+j
		print (c,j,dna[c])
for i in range (len(dna)):
	print (i, i % 3, dna[i])
"""
0 0 A
1 1 T
2 2 G
3 0 G
4 1 C
5 2 C
6 0 T
7 1 T
8 2 T
"""
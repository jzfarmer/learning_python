#!/usr/bin/env python3

import random
#random.seed(1) # comment-out this line to change sequence each time

# Write a program that stores random DNA sequence in a string
# The sequence should be 30 nt long
# On average, the sequence should be 60% AT
# Calculate the actual AT fraction while generating the sequence
# Report the length, AT fraction, and sequence

def at_seq(length, at):
	dna = ''
	for i in range(length):
	r = random.random()
	if r < t:
		r = random.random()
		if r < .5 : dna = 'A'
		else	  : dna = 'T'
	else:
		r = random.random()
		if r < .5: dna += 'G'
		else	 : dna += 'C'
	return dna

for i in range(3):
	print(at_seq(30, .6))
"""
python3 at_seq.py
30 0.6666666666666666 ATTACCGTAATCTACTATTAAGTCACAACC
"""

#!/usr/bin/env python3

# Write a program that computes the running sum from 1..n
# Also, have it compute the factorial of n while you're at it
# No, you may not use math.factorial()
# Use the same loop for both calculations

n = 5

def fac(n):
	fact = 1
	for i in range(1, n + 1):
		fact *= i
	return fact

for i in range(1,500):
	print(i, fac(i))
"""
python3 sumfac.py
5 15 120
"""

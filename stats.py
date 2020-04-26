#!/usr/bin/env python3

from math import sqrt
import fileinput

# Write a program that computes typical stats
# Count, Min, Max, Mean, Std. Dev, Median
# No, you cannot import any other modules!
data = []
for line in fileinput.input():
	if line.startswith('#'): continue
	line = line.rstrip() 
	data.append(float(line))
#median
data.sort()
if len(data) % 2 == 1: #odd
	median = index[len(data)/2]
else: #even
	m1 = int(len(data)/2)
	m2 = int(m1 - 1)
	median = (data[m1] + data[m2]) / 2	
#std dev
#sum
sum = 0
for x in data:
	sum += x
print(sum)
mean = sum/len(data)
variance = 0
for x in data:
	variance += (abs(mean - x))**2 / len(data)
stdev = sqrt(variance)
print(len(data), min(data), max(data), mean, f'{stdev:.3f}', median)
"""
python3 stats.py numbers.txt
Count: 10 
Minimum: -1.0
Maximum: 256.0
Mean: 29.147789999999997 #sum/total
Std. dev: 75.777 
Median 2.35914
"""

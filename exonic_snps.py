#!/usr/bin/env python3

# Write a program that intersects SNPs and coding exons
# snps.txt.gz has all the 23andme snps on chrom 21
# exons.txt.gz has all the coding exons on chrom 21
# Report:
#	number of SNPs
#	number of genes
#	number of exons
#	names of genes with SNPs in them
# Given that there are about 100x more genes in the human genome
#	And 100x more SNPs on the full chip
#	Estimate how long your program would take to run using all available data

import gzip 
import time
# get all the positions of the SNPs
t0 = time.perf_counter()
snps = []
with gzip.open('snps.txt.gz', 'rt') as fp:
	for line in fp.readlines():
		(id, chrom, pos, gtype) = line.split()
		snps.append(int(pos))
		
# read exons and xref to SNPs
genes = []
exons = 0
snp_genes = []
with gzip.open('exons.txt.gz', 'rt') as fp:
	for line in fp.readlines():
		(gene, beg, end) = line.split()
		exons += 1
		if gene not in genes: 
			genes.append(genes)
		beg = int(beg)
		end = int(end)
		print(beg)
		found = False
		for snp in snps:
			if snp >= beg and snp <= end:
				found = True
				break
		if found: 
			if gene not in snp_genes:
				snp_genes.append(gene)
snp_genes.sort()
t1 = time.perf_counter()
et = (t1 - t0) * 100 * 100 / 3600
#report section
print(f'SNPs:{len(snps)}')
print(f'Genes:{len(genes)}')
print(f'Exons:{exons}')
print(f'Genes w/ SNPs: {len(snp_genes)}')
print(f'Gene List: {snp_genes}')
print(f'Estimated Full Time: {et} hours')
"""
SNPs: 8607
Genes: 234
Exons: 2433
Genes w/ SNPs: 54
Gene List: ABCG1, AGPAT3, AIRE, AP000311.1, BACE2, BACH1, BRWD1, C21orf58, C2CD2, CBS, CHAF1B, CLDN17, COL6A1, COL6A2, DNMT3L, DOP1B, DSCAM, FAM243A, GART, IFNAR1, IFNGR2, ITGB2, KRTAP10-5, KRTAP10-7, KRTAP19-3, KRTAP25-1, MCM3AP, MORC3, PAXBP1, PCNT, PDE9A, PDXK, PIGP, PRDM15, PTTG1IP, PWP2, RRP1, RWDD2B, SCAF4, SIK1, SIM2, SLC37A1, SLC5A3, SOD1, SON, SYNJ1, TMPRSS15, TMPRSS2, TMPRSS3, TRAPPC10, TTC3, U2AF1, UMODL1, USP25
Estimated Full Time: 4.76 hours
"""

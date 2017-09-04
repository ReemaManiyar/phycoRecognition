### Import python modules and libraries
import pandas
import numpy
import glob
import os
import json
from Bio.SeqUtils import GC
from Bio import SeqIO
from Bio.SeqUtils import CodonUsage
#from __future__ import division
import math
import matplotlib
import re
from matplotlib import pyplot as plt
#import seaborn as sns
from scipy import stats
from scipy.stats.stats import pearsonr
### Load file containing list of plasmid and bacterial gene files
def codon_fre_extract():
	pos_neg_files = ['phycodnaviridae_virus_name.txt', 'phage_name.txt']
	data_path = 'C:\\Users\\Reema\\Documents\\SDSU_Education\\Thesis_Phyco\\'
	cdns = ['ATT', 'ATC', 'ATA', 'CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG', 'GTT', 'GTC', 'GTA', 'GTG', 'TTT', 'TTC', 'ATG', 'TGT', 'TGC', 'GCT', 'GCC', 'GCA', 'GCG', 'GGT', 'GGC', 'GGA', 'GGG', 'CCT', 'CCC', 'CCA', 'CCG', 'ACT', 'ACC', 'ACA', 'ACG', 'TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC', 'TAT', 'TAC', 'TGG', 'CAA', 'CAG', 'AAT', 'AAC', 'CAT', 'CAC', 'GAA', 'GAG', 'GAT', 'GAC', 'AAA', 'AAG', 'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG', 'TAA', 'TAG', 'TGA']
	codon_fre_dict = dict()
	for c in cdns:
		c_key = 'CDN_'+c
		codon_fre_dict[c_key] = []
	
	for example in pos_neg_files:
		example_path = data_path+example
		f = open(example_path, 'r')
		for dir in f.readlines():
			no_match = 0
			#print dir
			path= 'C:\\Users\\Reema\\Documents\\SDSU_Education\\Thesis_Phyco\\all_fna\\all.fna\\'
			path1= 'C:\\Users\\Reema\\Documents\\SDSU_Education\\Thesis_Phyco\\all_ffn\\all.ffn\\'
			dir = dir.strip('\n')
			dirpath= path+dir
			dirpath1= path1+dir
			#print "VIRUS = "+dirpath
			i= 'grinder-reads.fa'
			file_path= dirpath+'\\'+i
			file_handle = open(file_path, 'r')
			contig_num = 1
			for seq in SeqIO.parse(file_handle, "fasta"):
				num_gene = 0
				d= [0]*64
				inpstr=  seq.description
				complement_in_grinder = 'position=complement' in inpstr
				p = re.compile('[0-9]+\.\.[0-9]+')
				x = p.findall(inpstr)
				y= x[0].split('..')
				gmin = y[0]
				gmax = y[1]
				#print inpstr
				#print y
				matched = 0
				for j in os.listdir(dirpath1):
					if j.endswith(".ffn"):
						file_path1= dirpath1+'\\'+j
						file_handle1 = open(file_path1, 'r')
						for seq1 in SeqIO.parse(file_handle1, "fasta"):
							inpstr1= seq1.description
							q = re.compile('[0-9]+\-[0-9]+ ')
							z= q.findall(inpstr1)
							u= z[0].split('-')
							smin = u[0]
							smax = u[1]
							#print u
							complement_in_seq = '|:c' in inpstr1
							if  (complement_in_grinder and complement_in_seq) or (not(complement_in_grinder) and not (complement_in_seq)) :
								if gmin<= smin and gmax>=smax:
									na_seq = seq1.seq
									#GC_content = GC_content + GC(str(na_seq))
									#num_gene = num_gene +1
									#length= length + len(str(na_seq)
									codons = [str(na_seq)[n:n+3] for n in range(0,len(str(na_seq)),3)]
									codon_total = len(codons)
									c = list()
									for k in cdns:
										c.append(codons.count(k)/float(codon_total))
									d = numpy.add(d,c)
									matched = 1
									num_gene = num_gene +1
				if matched:
					d = d / float(num_gene)
					idx = 0
					#print d
					for c in cdns:
						c_key = 'CDN_'+c
						codon_fre_dict[c_key].append(d[idx])
						idx = idx + 1
				else:
					no_match = no_match + 1
					#print "Didnt match for ",contig_num
				contig_num = contig_num + 1
			
		f.close()	
	return codon_fre_dict

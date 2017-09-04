### Import python modules and libraries
import pandas
import numpy
import matplotlib
from matplotlib import pyplot as plt
import pylab
#import seaborn as sns
from Bio.SeqUtils import GC
from Bio import SeqIO
import os.path
import csv
import sys, getopt, csv
import re

### Load file listing plasmid and bacterial genomes. The file is a tab-separated file that contains the filename as the first column and the bacteria/plasmid name as the second columns
def gc_extract():
	pos_neg_files = ['phycodnaviridae_virus_name.txt', 'phage_name.txt']
	data_path = 'C:\\Users\\Reema\\Documents\\SDSU_Education\\Thesis_Phyco\\'
	gc_list = []
	num_gene_list = []
	gene_len_list = []
	name_list = []
	eg_list = ['pos', 'neg']
	eg_num = 0
	for example in pos_neg_files:
		prepend_str = eg_list[eg_num]
		example_path = data_path+example
		f = open(example_path, 'r')
		contig_num = 0
		for dir in f.readlines():
			#pvt_list = []
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
			for seq in SeqIO.parse(file_handle, "fasta"):
				GC_content = 0
				num_gene = 0
				length = 0
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
									matched = 1
									na_seq = seq1.seq
									GC_content = GC_content + GC(str(na_seq))
									num_gene = num_gene +1
									length= length + len(str(na_seq))
				if matched:
					GC_avg = round(GC_content/num_gene, 2)
					gc_list.append(GC_avg)
					avg_gene_length= length/num_gene
					gene_len_list.append(avg_gene_length)
					num_gene_list.append(num_gene)
					name_list.append(prepend_str+str(contig_num))
					contig_num = contig_num + 1
				else:
					no_match = no_match + 1
					#print "Didnt match for ",contig_num
			#print pvt_list
			#print "Not matched contigs = ", no_match
		f.close()
		eg_num = eg_num + 1
	ret_dict= {'Name': name_list, 'GCcontent': gc_list, 'num_of_gene': num_gene_list, 'length_of_gene': gene_len_list}
	return ret_dict	
	
			
			#	na_seq = seq.seq
			#	GC_content = GC_content + GC(str(na_seq))
			#	num_gene = num_gene +1
			#GC_avg = round(GC_content/num_gene, 2)
			#gc_list.append(GC_avg)

	




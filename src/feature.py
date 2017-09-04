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

##Import in-house modules
import GCcode
import codon_fre
import AA_fre
import con_seq
import pprint

#New empty dictionary
data_dict = dict()

#Call each data extractor here
#And assign it to dictionary

gc_dict = GCcode.gc_extract()
for k in gc_dict.keys():
	if not(k=='Name'):
		data_dict[k]= gc_dict[k]
name_list = gc_dict['Name']
#print name_list
print "GC is extracted..."

codon_fre_dict = codon_fre.codon_fre_extract()
for k in codon_fre_dict.keys():
	data_dict[k] = codon_fre_dict[k]
print "cdns freq are extracted..."

AAs_fre_dict = AA_fre.AA_fre_extract()
for k in AAs_fre_dict.keys():
	data_dict[k] = AAs_fre_dict[k]
print "AAs freq are extracted..."

#Conserved seq
print "Calling conserve sequence function..."
conservedseq = con_seq.con_seq()
for k in conservedseq.keys():
	data_dict[k] = conservedseq[k]
print "Conserved Seq are extracted..."

#Output to CSV
data_frame = pandas.DataFrame(data=data_dict, index=name_list)
data_frame.T.to_csv(path_or_buf='dataset.csv', sep=',')






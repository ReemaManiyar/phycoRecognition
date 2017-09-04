import random
### Import python modules and libraries
import pandas
import scipy
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix
import numpy
import matplotlib
from matplotlib import pyplot as plt
import pylab
#import seaborn as sns
from Bio.SeqUtils import GC
from Bio import SeqIO
import os.path
import csv
import pprint
from sklearn.externals import joblib

import sklearn.metrics
import matplotlib.pyplot as plt


path = "C:\\Users\\Reema\\Documents\\SDSU_Education\\Thesis_Phyco\\code\\dataset.csv"  

df = pandas.DataFrame.from_csv(path)
df1= df[:20]

"""
imp_aa= ['AAs_C','AAs_S','AAs_E','AAs_F','AAs_A','AAs_G','AAs_N']
df_impftr = pandas.DataFrame(index=index, columns=columns)
for i in range(0,len(imp_aa):
	df_impftr.loc[i]= df1.T[imp_aa[i]]	
"""	

df1.sort(ascending=True, inplace=True)
df2= df[20:84]
df2.sort(ascending=True, inplace=True)
df3= df[84:87]
df3.sort(ascending=True, inplace=True)
df4= df[87:]
df4.sort(ascending=True, inplace=True)
df_list =[df1,df2,df3,df4]
#df_list =[df1,df2,df3]

df =pandas.concat(df_list)

df = df.T
#print list(df.columns.values)
print "RFModel: Important DF shapre = ", df.shape
pos_eg_idx = [df.index.values[i] for i in range(0,len(df.index),1) if 'pos' in df.index.values[i]]
#neg_eg_idx = [df.index.values[i] for i in range(0,len(df.index),1) if 'neg' in df.index.values[i]]
total_pos =  len(pos_eg_idx)

pos_df = df[0:total_pos]
neg_df = df[total_pos:]
#print neg_df
#pos_train_df = pos_df.sample(n=10)
#print pos_train_df.index

pos_msk = numpy.random.rand(total_pos) < 0.7
pos_train_df = pos_df[pos_msk]
pos_test_df = pos_df[~pos_msk]
pos_train_target = [1]*len(pos_train_df)
pos_test_target = [1]*len(pos_test_df)

neg_msk = numpy.random.rand(len(neg_df)) < 0.7
neg_train_df = neg_df[neg_msk]
neg_test_df = neg_df[~neg_msk]
neg_train_target = [0]*len(neg_train_df)
neg_test_target = [0]*len(neg_test_df)

all_train_frame= [pos_train_df, neg_train_df]
all_train_df = pandas.concat(all_train_frame)

all_test_frame= [pos_test_df, neg_test_df]
all_test_df = pandas.concat(all_test_frame)

all_train_target = pos_train_target + neg_train_target
all_test_target = pos_test_target + neg_test_target
#print all_train_target
poslist =[]
acclist=[]
# random forest code
#for tree in range(20,200, 10): 
rf = RandomForestClassifier(n_estimators=150, min_samples_split=1, n_jobs=-1, oob_score=True, max_features="sqrt")
# fit the training data
#print('fitting the model')
rfObj = rf.fit(all_train_df, all_train_target)
pos_sc = rf.score(pos_test_df, pos_test_target)
all_sc = rf.score(all_test_df, all_test_target)
poslist.append(pos_sc)
acclist.append(all_sc)
	
print poslist
print acclist
	
#print "Number of features: ", rf.n_features_
	#oob_error = 1 - rf.oob_score_
fea_imp= rf.feature_importances_
#pprint.pprint(fea_imp)
	#print len(fea_imp)

#print rf.oob_score_
#print oob_error



# run model against test data
#predicted_probs = rf.predict_proba(all_test_df)
pos_sc = rf.score(pos_test_df, pos_test_target)
#falsepos_sc = rf.score(neg_test_df, pos_test_target)
pos_sc1 = rf.score(pos_train_df, pos_train_target)
neg_sc = rf.score(neg_test_df, neg_test_target)
neg_sc1 = rf.score(neg_train_df, neg_train_target)
#falseneg_sc = rf.score(pos_test_df, neg_test_target)
all_sc = rf.score(all_test_df, all_test_target)


print "================RESULT==================="
print "no of pos training samples: ", len(pos_train_df)
print "no of neg training samples: ", len(neg_train_df)
print "no of pos testing samples: ", len(pos_test_df)
print "no of neg traing samples: ",len(neg_test_df)
print "total testing samples: ", len(all_test_target)
print "total accuracy = ", all_sc
print "True Positive testing Rate = ",  pos_sc
print "True Positive training Rate = ",  pos_sc1
print "True Negative testing Rate = ", neg_sc
print "True Negative training Rate = ", neg_sc1
print "False Positive Rate = ", 1-neg_sc
print "False Negative Rate = ", 1-pos_sc

print "True Positive testing Rate = ",  pos_sc
#print "True Positive training Rate = ",  pos_sc1
print "True Negative Rate = ", neg_sc
### Calculations of Precision, Recall, F1, Accuracy
true_pos = pos_sc
false_pos= 1-neg_sc
false_neg= 1-pos_sc
prcsn= true_pos/(true_pos+false_pos)
recall= true_pos/(true_pos+false_neg)
f1_score= 2*prcsn*recall/(prcsn+recall)
print "================RESULT==================="
print "total accuracy = ", all_sc
print "precision = " ,  prcsn
print "recall = ", recall
print "F1 Score = ", f1_score

### Store fitted RF to a file --> phycoTool
joblib.dump(rf, 'phycoTool.pkl') 


#### ROC related
all_target = all_test_target + all_train_target
all_frame = [all_test_df, all_train_df]
all_df = pandas.concat(all_frame)
predProb = rfObj.predict_proba(all_test_df)
probs = numpy.amax(predProb, axis=1)
fpr, tpr, thresholds = sklearn.metrics.roc_curve(all_test_target , probs, pos_label=0)
roc_auc = sklearn.metrics.auc(fpr, tpr)
print "ROC area = ", roc_auc
plt.plot(fpr, tpr, lw=1, label='ROC')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic example')
plt.show()


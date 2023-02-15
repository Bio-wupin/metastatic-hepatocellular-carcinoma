#python
#this script is used to classify each event into "maintained, selected, unselected event"
#input data
#1.1 MT origin: paired PT~MT; paired RT~MT; paired MT~MT
##1.2 GISTIC result: sample level and tumor level
##1.3 Classification criteria: maintained, selected, unselected
#1.4 given a background null distribution (using passenger (non-driver) SCNAs)
#1.5 given driver SCNA eventa panel

import re
import os
import numpy as np

#Classification criteria: maintained, selected, unselected
print('read the Classification criteria!')
unselected=[['subclonal','absent'],['clonal','absent']]  #eg, subclonal in PT;  absent in MT
selected=[['absent','clonal'],['subclonal','clonal']]
maintained=[['subclonal','subclonal'],['clonal','clonal'],['clonal','subclonal']]

#GISTIC sample level result
print('read GISTIC result!')
dic_SCNA_sample={}#sample:[events]
path='/hpc/data/home/slst/zhangly/wupin/project/16-4th-sequencing-data-80011001-X3/6-GISTIC/3-GISTIC-output-center-default-median/all-sample-in-4th-EHM'
file=open(path+'/broad_values_by_arm.txt','r')
lines=file.readlines()
sample_index={}
data=lines[0].rstrip().split('\t')
for i in range(1,len(data)):
	sample_index.update({i:data[i]})#index:sample
	if data[i] not in dic_SCNA_sample.keys():
		dic_SCNA_sample.update({data[i]:[]})
for line in lines[1:]:
	data=line.rstrip().split('\t')
	for i in range(1,len(data)):
		sample=sample_index[i]
		if float(data[i])<=-0.5:
			SCNA=data[0]+'_Del'
			if SCNA not in dic_SCNA_sample[sample]:
				dic_SCNA_sample[sample].append(SCNA)
			else:
				print('error3')
		if float(data[i])>=0.5:
			SCNA=data[0]+'_Amp'
			if SCNA not in dic_SCNA_sample[sample]:
				dic_SCNA_sample[sample].append(SCNA)
			else:
				print('error4')

#GISTIC tumor level result
dic_SCNA_tumor={}#tumor:[events]
path='/hpc/data/home/slst/zhangly/wupin/project/16-4th-sequencing-data-80011001-X3/6-GISTIC/3-GISTIC-output-center-default-median/all-tumor-in-4th-EHM'
file=open(path+'/broad_values_by_arm.txt','r')
lines=file.readlines()
tumor_index={}
data=lines[0].rstrip().split('\t')
for i in range(1,len(data)):
	tumor_index.update({i:'_'.join(data[i].split('_')[0:2])})#index:tumor
	if '_'.join(data[i].split('_')[0:2]) not in dic_SCNA_tumor.keys():
		dic_SCNA_tumor.update({'_'.join(data[i].split('_')[0:2]):[]})
for line in lines[1:]:
	data=line.rstrip().split('\t')
	for i in range(1,len(data)):
		tumor=tumor_index[i]
		if float(data[i])<=-0.5:
			SCNA=data[0]+'_Del'
			if SCNA not in dic_SCNA_tumor[tumor]:
				dic_SCNA_tumor[tumor].append(SCNA)
			else:
				print('error1'+tumor+SCNA)
		if float(data[i])>=0.5:
			SCNA=data[0]+'_Amp'
			if SCNA not in dic_SCNA_tumor[tumor]:
				dic_SCNA_tumor[tumor].append(SCNA)
			else:
				print('error2'+tumor+SCNA)


#driver SCNA
driver_SCNA=[]
path='/hpc/data/home/slst/zhangly/wupin/project/binom-test-for-SCNA'
file=open(path+'/all-SCNA-in-HCC.txt','r')
lines=file.readlines()
for line in lines:
	data=line.rstrip()
	if data not in driver_SCNA:
		driver_SCNA.append(data)
print('all Arm-level SCNA in HCC: '+str(len(driver_SCNA)))


#MT origin: paired PT~MT; paired RT~MT; paired MT~MT
print('read all pairs: MT origin ~ MT')
all_pairs=[]
file=open('/hpc/data/home/slst/zhangly/wupin/project/binom-test-for-SCNA/MT-origin.txt','r')
lines=file.readlines()
for line in lines[1:]:
	data=line.rstrip().split('\t')
	tumor1=data[0]+'_'+data[1]
	tumor2=data[0]+'_'+data[2]
	if [tumor1,tumor2] not in all_pairs:
		all_pairs.append([tumor1,tumor2])

#calculate a background null distribution
print('calculate a background null distribution of proportions for all passenger SCNAs')
out=[]
for SCNA in driver_SCNA:
#for SCNA in driver_SCNA:
	N_main=0
	N_selected=0
	N_unselected=0
	for pair in all_pairs:
		tumor1=pair[0]
		tumor2=pair[1]
		tumor1_state=','
		tumor2_state=','
		#tumor1 state
		if SCNA in dic_SCNA_tumor[tumor1]:
			tumor1_state='clonal'
		else:
			#clonal/subclonal
			ll=[]
			for sample in dic_SCNA_sample.keys():
				if '_'.join(sample.split('_')[0:2])==tumor1:
					if SCNA in dic_SCNA_sample[sample]:
						ll.append('present')
					else:
						ll.append('absent')
			if len(ll)==0:
				print('error5')
			if ('present' in ll) and ('absent' in ll):
				tumor1_state='subclonal'
			if ('present' in ll) and ('absent' not in ll):
				tumor1_state='clonal'
			if ('present' not in ll) and ('absent' in ll):
				tumor1_state='absent'
		#tumor2 state
		if SCNA in dic_SCNA_tumor[tumor2]:
			tumor2_state='clonal'
		else:
			#clonal/subclonal
			ll=[]
			for sample in dic_SCNA_sample.keys():
				if '_'.join(sample.split('_')[0:2])==tumor2:
					if SCNA in dic_SCNA_sample[sample]:
						ll.append('present')
					else:
						ll.append('absent')
			if ('present' in ll) and ('absent' in ll):
				tumor2_state='subclonal'
			if ('present' in ll) and ('absent' not in ll):
				tumor2_state='clonal'
			if ('present' not in ll) and ('absent' in ll):
				tumor2_state='absent'
		if ',' in [tumor1_state,tumor2_state]:
			print('error6')
		#if ('clonal' in [tumor1_state,tumor2_state]) or ('subclonal' in [tumor1_state,tumor2_state]):
		#print(SCNA+'\t'+tumor1+'\t'+tumor1_state+'\t'+tumor2+'\t'+tumor2_state+'\n')
		if [tumor1_state,tumor2_state] in selected:
			#print(tumor1+'\t'+tumor1_state+'\t'+tumor2+'\t'+tumor2_state+'\t'+'selected')
			N_selected=N_selected+1
		if [tumor1_state,tumor2_state] in maintained:
			#print(tumor1+'\t'+tumor1_state+'\t'+tumor2+'\t'+tumor2_state+'\t'+'maintained')
			N_main=N_main+1
		if [tumor1_state,tumor2_state] in unselected:
			#print(tumor1+'\t'+tumor1_state+'\t'+tumor2+'\t'+tumor2_state+'\t'+'unselected')
			N_unselected=N_unselected+1
	print(SCNA+'\t'+'selected: '+str(N_selected)+'; unselected: '+str(N_unselected)+'; maintained: '+str(N_main))





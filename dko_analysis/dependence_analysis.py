#!/usr/bin/python

import os
import sys
import numpy




def compare_samples(sample):
	DIR='/mnt/thumper/home/snowdrop/TinyTree/Analysis/peak_find_macs2_low/'
	input1=open(DIR+sample+'_peaks.xls','r')
	all_input1=input1.readlines()
	output1=open(sample+'.bed','w')
	for i in all_input1[28:]:
		each=i.split()
		if float(each[6])>5:
			mid=(int(each[1])+int(each[2]))/2
			output1.write(each[0]+'\t'+str(max([0,mid-1500]))+'\t'+str(mid+1500)+'\t'+each[6]+'\n')
	output1.close()
	input1.close()

def fraction_h3k4me1_peaks():

	enhancer_dic={}
	input0=open('/mnt/thumper/home/snowdrop/TinyTree/Analysis/peak_define/clustering_heatmap_k4me1_k27ac/all_peak_den_merged.active.heatmap','r')
	all_input0=input0.readlines()
	for i in all_input0[1:]:
		each=i.split()
		enhancer_dic[each[0]]=1
	input0.close()
	input0=open('/mnt/thumper/home/snowdrop/TinyTree/Analysis/peak_define/clustering_heatmap_k4me1_k27ac/all_peak_den_merged.poised.heatmap','r')
	all_input0=input0.readlines()
	for i in all_input0[1:]:
		each=i.split()
		enhancer_dic[each[0]]=1
	input0.close()
	


	input1=open('../peak_define/all_peak_den_merged.txt','r')
	all_input1=input1.readlines()
	name=all_input1[0].split()
	i=0
	while i<len(name):
		if name[i]=='h3k4me1':
			h3k4me1_idx=i
		if name[i]=='h3k4me1_dko':
			h3k4me1_dko_idx=i
		if name[i]=='h3k4me3':
			h3k4me3_idx=i
		if name[i]=='h3k4me3_dko':
			h3k4me3_dko_idx=i
		i+=1
	
	count1=0
	count2=0
	for i in all_input1[1:]:
		each=i.split()
		id=each[0]+'.'+each[1]+'.'+each[2]
		### this is for the number of h3k4m3 decreased with decreased k4me1
		if float(each[h3k4me3_idx])>1:
			if float(each[h3k4me3_dko_idx])/float(each[h3k4me3_idx])<0.5:
				if float(each[h3k4me1_idx])>0.5:
					if float(each[h3k4me1_dko_idx])/float(each[h3k4me1_idx])<0.5:
						count1+=1
					else:
						count2+=1
				else:
					count2+=1
	print count1, count2
	
	count1=0
	count2=0
	count3=0
	count4=0
	for i in all_input1[1:]:
		each=i.split()
		id=each[0]+'.'+each[1]+'.'+each[2]
		### this is for the number of enhancer at h3k4me1 decreased
		if float(each[h3k4me1_idx])>1:
			if float(each[h3k4me1_dko_idx])/float(each[h3k4me1_idx])<0.5:
				if enhancer_dic.has_key(id):
					count1+=1
				else:
					count2+=1
			else:
				if enhancer_dic.has_key(id):
					count3+=1
				else:
					count4+=1
	print count1, count2, count3, count4


def fraction_class1_2_h3k4me1_peaks():

##class1: CR peaks overlapping with KMT2C/D independent h3k4me1 sites. class2: CR peaks overlapping with dependent h3k4me1 sites.
	input1=open('../peak_define/all_peak_den_merged.txt','r')
	all_input1=input1.readlines()
	name=all_input1[0].split()
	i=0
	while i<len(name):
		if name[i]=='h3k4me1':
			h3k4me1_idx=i
		if name[i]=='h3k4me1_dko':
			h3k4me1_dko_idx=i
		i+=1
	CR_list=[]
	for i in name[3:]:
		if len(i.split('_'))==1:
			CR_list.append(i)
	class1={}
	class2={}
	for i in all_input1[1:]:
		each=i.split()
		id=each[0]+'.'+each[1]+'.'+each[2]
		if float(each[h3k4me1_idx])>1:
			if float(each[h3k4me1_dko_idx])/float(each[h3k4me1_idx])<0.5:
				class2[id]=1
			else:
				class1[id]=1
	output1=open('class1_2_CR.txt','w')
	output1.write('CR\tclass1\tclass2\n')
	for i in CR_list:
		j=0
		while j<len(name):
			if name[j]==i:
				tmp_idx=j
			j+=1

		class1_count=0
		class2_count=0
		for j in all_input1[1:]:
			each=j.split()
			id=each[0]+'.'+each[1]+'.'+each[2]
			if float(each[tmp_idx])>numpy.log2(1.5):
				if class1.has_key(id):
					class1_count+=1
				elif class2.has_key(id):
					class2_count+=1
		output1.write(i+'\t'+str(class1_count)+'\t'+str(class2_count)+'\n')

def enhancer_associate_peak():
	input1=open('../peak_define/all_peak_den_merged.txt','r')
	all_input1=input1.readlines()
	output1=open('all_k4me1_k4me3_peak.bed','w')
	for i in all_input1[1:]:
		each=i.split()
		output1.write(each[0]+'\t'+each[1]+'\t'+each[2]+'\n')
	output1.close()
	os.system('~/software/BEDTools-Version-2.13.3/bin/intersectBed -a all_k4me1_k4me3_peak.bed -b /mnt/thumper/home/snowdrop/Data/Mouse/mm9/refseq.tss.25.25.bed -v > all_k4me1_k4me3_peak_enhancer.bed')

		

#fraction_class1_2_h3k4me1_peaks()
fraction_h3k4me1_peaks()
#enhancer_associate_peak()

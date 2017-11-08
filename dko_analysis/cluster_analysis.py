#!/usr/bin/python

import os
import numpy


def make_cluster():
	input1=open('../all_peak_den_merged.k4me1_k4me3.txt','r')
	all_input1=input1.readlines()
	output1=open('all_peak_den_merged.poised.heatmap','w')
	output2=open('all_peak_den_merged.active.heatmap','w')
	
	target_factor=['wstf','sur2','phf5a','med12','smc3','srsf1','srsf2','brg1','phrf1','baz1a','med1','snf2h','chd1']
	target_factor_idx=[]
	output1.write('id\t')
	output2.write('id\t')
	
	each=all_input1[0].split()
	for i in target_factor:
		j=0
		while j<len(each):
			if each[j]==i:
				target_factor_idx.append(j)
			j+=1
		
	for i in target_factor:
		output1.write(i+'\t')
		output2.write(i+'\t')
	output1.write('\n')
	output2.write('\n')
	for i in all_input1[1:]:
		each=i.split()
		if (float(each[3])>=1 and float(each[4])<0.263):
			if float(each[10])>=1:
				id=each[0]+'.'+each[1]+'.'+each[2]
				output2.write(id+'\t')
				for j in target_factor_idx:
					output2.write(each[j]+'\t')
				output2.write('\n')
			else:
				id=each[0]+'.'+each[1]+'.'+each[2]
				output1.write(id+'\t')
				for j in target_factor_idx:
					output1.write(each[j]+'\t')
				output1.write('\n')
	

	input1.close()
	output1.close()
	output2.close()

def make_clustering():
	os.system('~/software/cluster-1.52a/bin/cluster -f all_peak_den_merged.poised.heatmap -g 7 -k 15')
	os.system('~/software/cluster-1.52a/bin/cluster -f all_peak_den_merged.active.heatmap -g 7 -k 15')

def change_format_qq_data():
	input1=open('qq_CR_data.heatmap','r')
	input2=open('all_peak_den_merged.k4me1_k4me3.heatmap','r')
	output1=open('qq_CR_data.heatmap.label','w')
	all_input1=input1.readlines()
	all_input2=input2.readlines()
	output1.write('id\t')
	for i in all_input2[0].split()[3:]:
		output1.write(i+'\t')
	output1.write('\n')
	i=1
	while i<len(all_input1):
		id=all_input2[i].split()[0]
		output1.write(id+'\t')
		value=all_input1[i].split()
		for j in value[1:]:
			output1.write(j+'\t')
		output1.write('\n')
		i+=1

def count_cluster():
	input1=open('qq_CR_data.heatmap_K_G30.kgg','r')
	all_input1=input1.readlines()
	all_count={}
	for i in all_input1[1:]:
		each=i.split()
		if all_count.has_key(each[1]):
			all_count[each[1]]+=1
		else:
			all_count[each[1]]=1
	for i in range(30):
		print i,all_count[str(i)]

def reorder_cluster():
	cluster_order=[4,6,11,1,12,7,3,13,14,10,9,2,8,5,0]
	input1=open('all_peak_den_merged.active_K_G15.cdt','r')
	
	#cluster_order=[0,2,3,7,9,5,6,1,8,13,12,11,10,4,14]
	#input1=open('all_peak_den_merged.poised_K_G15.cdt','r')
	
	all_input1=input1.readlines()
	all_data={}
	for i in all_input1[1:]:
		each=i.split()
		all_data[each[0]]=i
	input2=open('all_peak_den_merged.active_K_G15.kgg','r')
	all_input2=input2.readlines()
	output1=open('all_peak_den_merged.active.order.cdt','w')
	
	#input2=open('all_peak_den_merged.poised_K_G15.kgg','r')
	#all_input2=input2.readlines()
	#output1=open('all_peak_den_merged.poised.order.cdt','w')
	
	output1.write('id\t')
	i=1
	header=all_input1[0].split()
	while i<len(header):
		output1.write(header[i]+'\t')
		i+=1
	output1.write('\n')
	
	for i in cluster_order:
		for j in all_input2[1:]:
			each=j.split()
			if each[1]==str(i):
				output1.write(all_data[each[0]])


#make_cluster()
#make_clustering()
#change_format_qq_data()
#count_cluster()
#reorder_cluster()

#!/usr/bin/python

import numpy
import operator

def define_CR_den_k4me1_k4me3():
# define peak density by taking the average of any replicates with over 1.5 fold enrichment over input
# otherwise define as zero
	input1=open('RPKM/brg1.rep1.filter.bam.h3k4me1.3.method1.peak.rpkm','r')
	all_input1=input1.readlines()
	peak_dic={}
	for i in all_input1:
		each=i.split()
		id=each[0]+'.'+each[1]+'.'+each[2]
		if i.count('h3k4me3')>0:
			if i.count('h3k4me1')>0:
				peak_dic[id]='h3k4me1_3'
			else:
				peak_dic[id]='h3k4me3'
		else:
			peak_dic[id]='h3k4me1'
	input1.close()

	output1=open('all_peak_den_merged.txt','w')
	output1.write('chr\tstart\tend\tpeak_type\t')
	sample_list=['baf155.rep1.filter.bam','baf155.rep2.filter.bam']
        current_sample='XXX'
	input1=open('/mnt/thumper/home/snowdrop/TinyTree/Data/all_filter.bam.txt.gb','r')  # select samples after checking manually
        all_input1=input1.readlines()
	#all_input1.append('test')
        all_peak_den=[]
	for i in all_input1:
        	each=i.split()
        	sample=each[0].split('.')[0]
        	if sample!=current_sample:
			print sample
			sample_peak_den=[]
			input2=open('RPKM/'+sample+'.normal.rpkm','r')
			output1.write(sample+'\t')
			current_sample=sample
			all_input2=input2.readlines()
			total_sample=len(all_input2[0].split())-3
			max_sample_idx=3
			max_peak=0
			j=3
			while j<len(all_input2[0].split()):
				count_peak=0
				for k in all_input2[1:]:
					each=k.split()
					if float(each[j])>numpy.log2(1.5):
						count_peak+=1	
					if count_peak>max_peak:
						max_peak=count_peak
						max_sample_idx=j	
				j+=1

			for j in all_input2[1:]:
				each=j.split()
				if len(each)>4:
					if float(each[max_sample_idx])>numpy.log2(1.5):
						peak_den=[]
						for k in each[3:]:
							if float(k)>numpy.log2(1.5):
								peak_den.append(float(k))
						if len(peak_den)>1:
							peak_den=sum(peak_den)/float(len(peak_den))
						else:
							peak_den=0
					else:
						peak_den=0
				else:
					if float(each[max_sample_idx])>numpy.log2(1.5):
						peak_den=each[max_sample_idx]
					else:
						peak_den=0	
				sample_peak_den.append(peak_den)	
			all_peak_den.append(sample_peak_den)	
	output1.write('\n')
	i=1
	while i<len(all_input2):
		each=all_input2[i].split()
		id=each[0]+'.'+each[1]+'.'+each[2]
		type=peak_dic[id]
		output1.write(each[0]+'\t'+each[1]+'\t'+each[2]+'\t'+type+'\t')
		#print all_peak_den
		for j in all_peak_den:
			output1.write(str(j[i-1])+'\t')
		output1.write('\n')
		i+=1
	


#define_CR_den_k4me1_k4me3()

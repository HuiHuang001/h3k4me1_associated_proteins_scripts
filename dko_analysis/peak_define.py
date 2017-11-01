#!/usr/bin/python

import operator
import os
import commands
import fileinput
import numpy

from functools import partial
from multiprocessing.dummy import Pool
from subprocess import call



def callpeak():


        commands=[]

        input1=open('/mnt/silencer2/home/huh025/AL_HH/merged.bam.txt', 'r')
        ## change bam file list to call peaks for merged sample as well as individual replicates 
        all_input=input1.readlines()
        
        for i in all_input:
                i=i.strip('\n')
                each=i.split('nput.')
                
                commands.append('macs2 callpeak -c /mnt/silencer2/home/huh025/AL_HH/merged_bam/'+i+' -t /mnt/silencer2/home/huh025/AL_HH/merged_bam/'+each[1]+' -f BAM -g mm -n '+each[1]+' -B -p 1e-5 -m 5 50 --outdir /mnt/silencer2/home/huh025/AL_HH/merged_peak')
               
        print commands       
        pool = Pool(15) # two concurrent commands at a time
        for i, returncode in enumerate(pool.imap(partial(call, shell=True), commands)):
               if returncode != 0:
                       	print("%d command failed: %d" % (i, returncode))
                	
                        
##callpeak()



def run_program():
	#sample_list=['baf155.rep1','baf155.rep2','baf155.rep3','baf155.rep4']
	sample_list=['baf155.rep3']
	current_sample='baf155'
	
	input1=open('/mnt/silencer2/home/huh025/AL_HH/good_peak.list.txt','r')  # sort list by sample name
   	all_input1=input1.readlines()
   	for i in all_input1[1:]:
    		i=i.strip('\n')
		each=i.split()
		sample=each[0].split('.')[0]
		if sample!=current_sample:
		    		sample_list.append(current_sample+'.merged.bam')
		#define peaks for current_sample
				if len(sample_list)==2:
						print sample_list[0]
						os.system('cp /mnt/silencer2/home/huh025/AL_HH/peak/'+sample_list[0]+'_peaks.narrowPeak '+'/mnt/silencer2/home/huh025/AL_HH/method1_peak/'+sample_list[0]+'.method1.peak')
						sample_list=[]
                        current_sample=sample
                        sample_list.append(i.split('_peaks')[0])
				else:
						method1(sample_list)			
						sample_list=[]
						current_sample=sample
						sample_list.append(i.split('_peaks')[0])
		else:
		#load replicates in sample_list
			sample_list.append(i.split('_peaks')[0])
	sample_list.append(current_sample+'.merged.bam')
	method1(sample_list)
	

def method1(sample_list):
## Due to limited quality of some TF ChIP, I define peaks by pooled replicates, keep the peaks that are reproducible in at least one replicate.
	
	
	
	        save_file=sample_list[0].split('.')[0]+'.method1.peak'
	        merge_file=sample_list[-1]
	        print merge_file
	        for i in sample_list:
	        #define reproducible peaks of the most peaked replicate
		        if i!=merge_file:
			    		os.system('bedtools intersect -a /mnt/silencer2/home/huh025/AL_HH/merged_peak/'+merge_file+'_peaks.narrowPeak'+' -b /mnt/silencer2/home/huh025/AL_HH/peak/'+i+'_peaks.narrowPeak -u >> tmp.txt')
		os.system('sort -k1,1 -k2,2n tmp.txt | uniq > '+'/mnt/silencer2/home/huh025/AL_HH/method1_peak/'+save_file)
	        os.system('rm tmp.txt')


def check_peak_num():
	input1=open('list.txt','r')
	all_input1=input1.readlines()
	for i in all_input1:
		each=i.split()
		status,output=commands.getstatusoutput('cat '+each[0]+' | wc -l')
		print each[0], output

def merge_h3k4me1_2_3_peak():
##merge reproducible h3k4me1, h3k4me2, h3k4me3 peaks.
	os.system('cat /mnt/silencer2/home/huh025/AL_HH/method1_peak/h3k4me1.method1.peak /mnt/silencer2/home/huh025/AL_HH/method1_peak/h3k4me2.method1.peak /mnt/silencer2/home/huh025/AL_HH/method1_peak/h3k4me3.method1.peak | sort -k1,1 -k2,2n - > tmp.bed')
	os.system('bedtools merge -d 1000 -c 4 -o collapse -i tmp.bed > h3k4me1.2.3.method1.peak')



#merge_h3k4me1_3_peak()






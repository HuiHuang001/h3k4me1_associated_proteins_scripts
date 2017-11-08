#!/usr/bin/pyhon
#due to the challenge of doing chip-seq with TFs, use method1 to find confident peaks.
# method1: find the best rep in terms of number of peak and then check reproducible peak from other replicates


import operator
import os
import commands
import fileinput
import numpy

def run_program():
	#sample_list=['baf155.rep1','baf155.rep2','baf155.rep3','baf155.rep4']
	sample_list=['baf155.rep1','baf155.rep2']
	current_sample='baf155'
	#input1=open('/mnt/thumper/home/snowdrop/TinyTree/Data/all_filter.bam.txt','r')
	input1=open('/mnt/thumper/home/snowdrop/TinyTree/Data/all_filter.bam.txt.gb','r')  # select samples after checking manually
        all_input1=input1.readlines()
        for i in all_input1[2:]:
		each=i.split()
		sample=each[0].split('.')[0]
		if sample!=current_sample:
			#print sample_list
			if len(sample_list)==1:
				print sample_list[0]
				os.system('cp ../peak_find_macs2/'+sample_list[0]+'_peaks.narrowPeak '+sample_list[0]+'.method1.peak')
			else:
				method1(sample_list)			
			sample_list=[]
			current_sample=sample
			sample_list.append(i.split('.filter')[0])
		else:
			sample_list.append(i.split('.filter')[0])

	#print sample_list
	if len(sample_list)==1:
		os.system('cp ../peak_find_macs2/'+sample_list[0]+'_peaks.narrowPeak '+sample_list[0]+'.method1.peak')
	else:
		method1(sample_list)

def method1(sample_list):
	#print sample_list
	# find the most peaked sample
	max_file_name=''
	max_peak_num=0
	for i in sample_list:
		status,output=commands.getstatusoutput('cat ../peak_find_macs2/'+i+'_peaks.narrowPeak | wc -l')
		if int(output)>max_peak_num:
			max_peak_num=int(output)
			max_file_name=i
			save_file=max_file_name.split('.')[0]+'.method1.peak'
	
	print max_file_name
	#for i in sample_list:
	#	if i!=max_file_name:
	#		os.system('/mnt/thumper/home/snowdrop/software/BEDTools-Version-2.13.3/bin/intersectBed -a ../peak_find_macs2/'+max_file_name+'_peaks.narrowPeak -b ../peak_find_macs2_low/'+i+'_peaks.narrowPeak -u >> tmp.txt')

	#os.system('sort -k1,1 -k2,2n tmp.txt | uniq > '+save_file)		
	#os.system('rm tmp.txt')


def check_peak_num():
	input1=open('list.txt','r')
	all_input1=input1.readlines()
	for i in all_input1:
		each=i.split()
		status,output=commands.getstatusoutput('cat '+each[0]+' | wc -l')
		print each[0], output

def merge_h3k4me1_3_peak():
	os.system('cat h3k4me1.method1.peak h3k4me3.method1.peak | sort -k1,1 -k2,2n - > tmp.bed')
	os.system('/mnt/thumper/home/snowdrop/software/BEDTools-Version-2.13.3/bin/mergeBed -d 1000 -nms -i tmp.bed > h3k4me1.3.method1.peak')


def rpkm_values():
	input1=open('/mnt/thumper/home/snowdrop/TinyTree/Data/all_filter.bam.txt.gb','r')  # select samples after checking manually
        all_input1=input1.readlines()
	for i in all_input1[124:125]:
		each=i.split()
		print each[0]
		os.system('/mnt/thumper/home/snowdrop/software/BEDTools-Version-2.13.3/bin/coverageBed -abam /mnt/thumper/home/snowdrop/TinyTree/Data/'+each[0]+' -b h3k4me1.3.method1.peak > RPKM/'+each[0]+'.h3k4me1.3.method1.peak.rpkm')
	
	#os.system('/mnt/thumper/home/snowdrop/software/BEDTools-Version-2.13.3/bin/coverageBed -abam /mnt/thumper/home/snowdrop/TinyTree/Data/merge_data/input_wt.all.bam -b h3k4me1.3.method1.peak > RPKM/input_wt.all.h3k4me1.3.method1.peak.rpkm')
	#os.system('/mnt/thumper/home/snowdrop/software/BEDTools-Version-2.13.3/bin/coverageBed -abam /mnt/thumper/home/snowdrop/TinyTree/Data/merge_data/input_dko.all.bam -b h3k4me1.3.method1.peak > RPKM/input_dko.all.h3k4me1.3.method1.peak.rpkm')
	#os.system('/mnt/thumper/home/snowdrop/software/BEDTools-Version-2.13.3/bin/coverageBed -abam /mnt/thumper/home/snowdrop/TinyTree/Data/merge_data/input.all.bam -b h3k4me1.3.method1.peak > RPKM/input.all.h3k4me1.3.method1.peak.rpkm')

def line_num_bam_file():
	input1=open('/mnt/thumper/home/snowdrop/TinyTree/Data/all_filter.bam.txt.gb','r')  # select samples after checking manually
        all_input1=input1.readlines()
	for i in all_input1[124:]:
		each=i.split()
		status,output=commands.getstatusoutput('samtools view /mnt/thumper/home/snowdrop/TinyTree/Data/'+each[0]+' | wc -l')
		print each[0],output
	#for i in ['input.all.bam','input_wt.all.bam','input_dko.all.bam']:	
	#	status,output=commands.getstatusoutput('samtools view /mnt/thumper/home/snowdrop/TinyTree/Data/merge_data/'+i+' | wc -l')
	#	print each[0],output
	
def make_rpkm(file_list):
	input1=open('bam_file_read_number.txt','r')
	read_dic={}
	all_input1=input1.readlines()
	for i in all_input1:
		each=i.split()
		read_dic[each[0]]=int(each[1])
		# snf2h.rep1.filter.bam 30888773
	#[baf155.rep1,baf155,rep2,input.all.bam]
	input1.close()

	print file_list
	input1=open('RPKM/'+file_list[-1]+'.h3k4me1.3.method1.peak.rpkm','r')
	all_input1=input1.readlines()
	input_rpkm=[]
	total_input=read_dic[file_list[-1]]
	for i in all_input1:
		each=i.split()
		dist=float(each[6])
                tag=float(each[4])
                rpkm=tag/(dist/1000.0)/(total_input/1000000.0)
                input_rpkm.append([each[0],each[1],each[2],rpkm])
	input1.close()
	all_sample_rpkm=[]
	for i in file_list[:-1]:
		sample=i
		total_sample=read_dic[sample]
		input1=open('RPKM/'+sample+'.h3k4me1.3.method1.peak.rpkm','r')
		all_input1=input1.readlines()
		sample_rpkm=[]
		for j in all_input1:
			each=j.split()
			dist=float(each[6])
                	tag=float(each[4])
                	rpkm=tag/(dist/1000.0)/(total_sample/1000000.0)
                	sample_rpkm.append([each[0],each[1],each[2],rpkm])
		input1.close()
		all_sample_rpkm.append(sample_rpkm)
	
	output1=open('RPKM/'+file_list[0].split('.')[0]+'.normal.rpkm','w')
	output1.write('chr\tstart\tend\t')
	for i in file_list[:-1]:
		output1.write(i+'\t')
	output1.write('\n')
	i=0
	while i<len(input_rpkm):
		output1.write(input_rpkm[i][0]+'\t'+str(input_rpkm[i][1])+'\t'+str(input_rpkm[i][2])+'\t')
		for j in all_sample_rpkm:
			output1.write(str(numpy.log2((j[i][3]+0.1)/(input_rpkm[i][3]+0.1)))+'\t')
		output1.write('\n')
		i+=1
	output1.close()
		
		
	
def rpkm_values_all_sample():
	sample_list=['baf155.rep1.filter.bam','baf155.rep2.filter.bam']
	current_sample='baf155'
	input1=open('/mnt/thumper/home/snowdrop/TinyTree/Data/all_filter.bam.txt.gb','r')  # select samples after checking manually
        all_input1=input1.readlines()
        for i in all_input1[124:]:
		each=i.split()
		sample=each[0].split('.')[0]
		if sample!=current_sample:
			if current_sample.count('wt')==1:
				sample_list.append('input_wt.all')
			elif current_sample.count('dko')==1:			
				sample_list.append('input_dko.all')
			else:
				sample_list.append('input.all')

			make_rpkm(sample_list)
			sample_list=[]
			current_sample=sample
			sample_list.append(each[0])
		else:
			sample_list.append(each[0])
		
	make_rpkm(sample_list)

def merge_sample():
	sample_list=['baf155.rep1.filter.bam','baf155.rep2.filter.bam']
	current_sample='baf155'
	input1=open('/mnt/thumper/home/snowdrop/TinyTree/Data/all_filter.bam.txt.gb','r')  # select samples after checking manually
        all_input1=input1.readlines()
        for i in all_input1[2:]:
		each=i.split()
		sample=each[0].split('.')[0]
		if sample!=current_sample:

			if len(sample_list)==1:
				script='cp ../'+sample_list[0]+' '+current_sample+'.all.bam'
				print script
			else:	

				script='samtools merge '+current_sample+'.all.bam'
				for j in sample_list:
					script=script+' ../'+j

			#print script+' &'

			sample_list=[]
			current_sample=sample
			sample_list.append(each[0])
		else:
			sample_list.append(each[0])
	print sample_list	





#run_program()
#check_peak_num()
#rpkm_values()
#line_num_bam_file()
#rpkm_values_all_sample()


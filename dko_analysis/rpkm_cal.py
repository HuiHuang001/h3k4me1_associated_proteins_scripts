#!/usr/bin/python

import operator
import os
import commands
import fileinput
import numpy

def sample_coverage():

#count coverage for each sample at h3k4me1,h3k4me2, h3k4me3 merged peaks.
	input1=open('/mnt/silencer2/home/huh025/AL_HH/input.peak.txt','r')  	# calculate coverage for all samples including input
        all_input1=input1.readlines()
        commands=[]
		
	for i in all_input1:
		        each=i.split('_peaks')
		        print each[0]
		        command='bedtools coverage -b /mnt/silencer2/home/huh025/AL_HH/bam/'+each[0]+'.filt.nodup.srt.bam'+' -a /mnt/silencer2/home/huh025/AL_HH/h3k4me1.2.3.peak.1kb.bin.bed > /mnt/silencer2/home/huh025/AL_HH/RPKM_mm9_10kb_bin/'+each[0]+'h3k4me1.2.3.peak.1kb.bin.bed.rpkm'
		        print command
			commands.append(command)
	
        pool = Pool(8) # 8 concurrent commands at a time
        for i, returncode in enumerate(pool.imap(partial(call, shell=True), commands)):
        	if returncode != 0:
        		print("%d command failed: %d" % (i, returncode))
	


def make_rpkm(file_list):
#input read number file
	input1=open('/mnt/silencer2/home/huh025/AL_HH/bam.read.num.txt','r')
	read_dic={}
	all_input1=input1.readlines()
	for i in all_input1:
		each=i.split()
		read_dic[each[0]]=int(each[1])
		# snf2h.rep1.filter.bam 30888773
	#[baf155.rep1,baf155,rep2,input.all.bam]
	input1.close()

	print file_list
	#make a dictrary of sample and input 
	input1=open('/mnt/silencer2/home/huh025/AL_HH/chip.info.txt','r')
	sample_dic={}
	all_input1=input1.readlines()
	for i in all_input1:
		each=i.split()
		sample_dic[each[0]]=each[1]
		# snf2h.rep1.filter.bam 30888773
	#[baf155.rep1,baf155,rep2,input.all.bam]
	input1.close()
	
	
	all_sample_rpkm=[]
	all_input_rpkm=[]
	for i in file_list:
	
		sample=i
		input_file=i.split('.')
#open coverage file of input	    
	    	input2=open('/mnt/silencer2/home/huh025/AL_HH/RPKM_123/'+sample_dic[i].split('.filt')[0]+'.h3k4me1.2.3.method1.peak.rpkm','r')
	    	all_input2=input2.readlines()
	    	input_rpkm=[]
	    	total_input=read_dic[sample_dic[i]]
	    	for j in all_input2:
			each=j.split()
			#print each
			dist=float(each[6])
            		tag=float(each[4])
            		rpkm=tag/(dist/1000.0)/(total_input/1000000.0)
            		input_rpkm.append([each[0],each[1],each[2],rpkm])
		input2.close()
		all_input_rpkm.append(input_rpkm)
	
	
		total_sample=read_dic[sample]
		input3=open('/mnt/silencer2/home/huh025/AL_HH/RPKM_123/'+sample.split('.filt')[0]+'.h3k4me1.2.3.method1.peak.rpkm','r')
		all_input3=input3.readlines()
		sample_rpkm=[]
		for j in all_input3:
			
			each=j.split()
			dist=float(each[6])
                	tag=float(each[4])
               		rpkm=tag/(dist/1000.0)/(total_sample/1000000.0)
                	sample_rpkm.append([each[0],each[1],each[2],rpkm])
		input3.close()
		all_sample_rpkm.append(sample_rpkm)
	output1=open('/mnt/silencer2/home/huh025/AL_HH/tmp123.txt','w')
	i=0
	while i<len(input_rpkm):
		output1.write(input_rpkm[i][0]+'\t'+str(input_rpkm[i][1])+'\t'+str(input_rpkm[i][2])+'\t')
		for j in all_sample_rpkm:
		#write line i of all_sample to output line i
			#print j[i]
			output1.write(str(j[i][3]+0.1)+'\t')

		for k in all_input_rpkm:
			#print k[i]

		    	output1.write(str(k[i][3]+0.1)+'\t')
		output1.write('\n')
		i+=1
	output1.close()
	l=len(file_list)
	input4=open('/mnt/silencer2/home/huh025/AL_HH/tmp123.txt','r')
	output2=open('/mnt/silencer2/home/huh025/AL_HH/RPKM_123/'+file_list[0].split('.')[0]+'.normal.rpkm','w')
	
	#write header chr start end sample
	output2.write('chr\tstart\tend\t')
	for i in file_list:
		output2.write(i.split('.')[0]+i.split('.')[1]+'\t')
	output2.write('\n')
	
	i=0
	all_input4=input4.readlines()
	#while i<len(input_rpkm):
		#output2.write(input_rpkm[i][0]+'\t'+str(input_rpkm[i][1])+'\t'+str(input_rpkm[i][2])+'\t')
	
	for j in all_input4:
	        
	        	each=j.split()
	       		#print each
			output2.write(each[0]+'\t'+each[1]+'\t'+each[2]+'\t')
			i=0
			while i<len(file_list): 
	        		output2.write(str(numpy.log2(float(each[3+i])/float(each[3+i+l])))+'\t')
				i=i+1
	    		output2.write('\n')
		
		
	
def rpkm_values_all_sample():
	sample_list=['baf155.rep3.filt.nodup.srt.bam']
	current_sample='baf155'
	input1=open('/mnt/silencer2/home/huh025/AL_HH/chip.info.txt','r')  # select samples after checking manually
        all_input1=input1.readlines()
	for i in all_input1[1:]:
		        each=i.split()
		        sample=each[0].split('.')[0]
			if sample!=current_sample:
		             #   sample_list.append(current_sample+'.all.bam')
			    	#if current_sample.count('wt')==1:
				   #     sample_list.append('input_wt.all.bam')
			    	#elif current_sample.count('dko')==1:							
					#sample_list.append('input_dko.all.bam')
			    	#else:
				#sample_list.append('input.all.bam')

			        make_rpkm(sample_list)
		                sample_list=[]
			        current_sample=sample
			        sample_list.append(each[0])
			else:
			            sample_list.append(each[0])
	#sample_list.append(current_sample+'.all.bam')
	#if current_sample.count('wt')==1:
	#	sample_list.append('input_wt.all.bam')
    #    elif current_sample.count('dko')==1:                                                    
     #   	sample_list.append('input_dko.all.bam')
     #   else:
		#sample_list.append('input.all.bam')
	make_rpkm(sample_list)
	
	
	

rpkm_values_all_sample()





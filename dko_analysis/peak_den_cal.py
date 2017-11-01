
#!/usr/bin/python

import numpy
import operator

def define_CR_den_k4me1_k4me3():
# define peak density by taking the average of any sample over 1.5-fold enrichment for CR, 2 fold for H3K27ac 
# otherwise define as zero
	input1=open('/mnt/silencer2/home/huh025/AL_HH/RPKM_123/baf155.rep3.h3k4me1.2.3.method1.peak.rpkm','r')
	all_input1=input1.readlines()
	peak_dic={}
	for i in all_input1:
		each=i.split()
		id=each[0]+'.'+each[1]+'.'+each[2]
		if i.count('h3k4me1')>0:
			if i.count('h3k4me2')>0:
				if i.count('h3k4me3')>0:
					peak_dic[id]='h3k4me1_2_3'
				else:
					peak_dic[id]='h3k4me1_2'
			else:
				if i.count('h3k4me3')>0:
					peak_dic[id]='h3k4me1_3'
				else:
					peak_dic[id]='h3k4me1'
		else:
			if i.count('h3k4me2')>0:
                                if i.count('h3k4me3')>0:
                                        peak_dic[id]='h3k4me2_3'
                                else:
                                        peak_dic[id]='h3k4me2'
                        else:
                                if i.count('h3k4me3')>0:
                                        peak_dic[id]='h3k4me3'
	input1.close()

	output1=open('/mnt/silencer2/home/huh025/AL_HH/all_peak_den_123_merged.txt','w')
	#average log2 fold enrichment over input of best replicates for every sample 
	output1.write('chr\tstart\tend\tpeak_type\t')
	#sample_list=['baf155.rep1.filt.nodup.srt.bam','baf155.rep2.filt.nodup.srt.bam']
    	#current_sample='baf155'
	input1=open('/mnt/silencer2/home/huh025/AL_HH/all.normal.list.txt','r')  # select samples after checking manually
    	all_input1=input1.readlines()
	
    	all_peak_den=[]
	for i in all_input1:
        		
			each=i.split()
        		sample=each[0].split('.')[0]
        		
			print sample
			sample_peak_den=[]
			input2=open('/mnt/silencer2/home/huh025/AL_HH/RPKM_123/'+sample+'.normal.rpkm','r')
			#log2 fold enrichment over merged input for all the replicates, including merged file. 
			output1.write(sample+'\t')
			
			all_input2=input2.readlines()
			total_sample=len(all_input2[0].split())-3
			
		
			for j in all_input2[1:]:
				each=j.split()
				
				if len(each)>4:
				#only calculate peak_den when it contains more than two replicates
					
						peak_den=[]

						for k in each[3:]:
							if float(k)>numpy.log2(1.5):
								peak_den.append(float(k))
							#else:
							#	peak_den=0
						if len(peak_den)>1:
							peak_den=sum(peak_den)/float(len(peak_den))
				            		#only calculate peak_den when more than two replicates have >1.5 fold enrichment.
						else:
							peak_den=0
				
				else:
						peak_den=[]
						if float(each[3])>numpy.log(1.5):
							peak_den=float(each[3])
						else:
							peak_den=0
				sample_peak_den.append(peak_den)	
				#inner for loop record peak density for current sample in sample_peak_den
			all_peak_den.append(sample_peak_den)
			#outer for loop record all sample peak density in all_peak_den
				
	output1.write('\n')
	i=1
	while i<len(all_input2):
	#for all lines
		each=all_input2[i].split()
		id=each[0]+'.'+each[1]+'.'+each[2]
		type=peak_dic[id]
		output1.write(each[0]+'\t'+each[1]+'\t'+each[2]+'\t'+type+'\t')
		#print all_peak_den
		for j in all_peak_den:
		#print line i for every sample
			output1.write(str(j[i-1])+'\t')
		output1.write('\n')
		i+=1
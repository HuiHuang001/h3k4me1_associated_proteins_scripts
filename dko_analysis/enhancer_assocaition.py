#!/usr/bin/python

import operator
import os
import commands
import fileinput
import numpy

def distal_h3k4me1_2_3_peak():
# define regions +-2.5 kb away from tss
		input1=open('/mnt/silencer2/home/huh025/AL_HH/all_peak_den_123_merged.txt','r')
		all_input1=input1.readlines()
		output1=open('/mnt/silencer2/home/huh025/AL_HH/h3k4me1_2_3_peak.bed','w')
		#find h3k4me1 peaks
		for i in all_input1[1:]:
			each=i.split()
			if each[3].count('1')>0:
				output1.write(each[0]+'\t'+each[1]+'\t'+each[2]+'\n')
		output1.close()
		os.system('bedtools intersect -a /mnt/silencer2/home/huh025/AL_HH/h3k4me1_2_3_peak.bed -b /mnt/silencer2/home/huh025/AL_HH/mm9.gencode.25tss25.bed -v > /mnt/silencer2/home/huh025/AL_HH/distal_h3k4me1_2_3_peak.bed')
		input2=open('/mnt/silencer2/home/huh025/AL_HH/distal_h3k4me1_2_3_peak.bed','r')
		all_input2=input2.readlines()
		distal_dic={}
		##define distal h3k4me1_2_3 peaks overlapping with h3k27ac as active enhancer, non-overlapping ones as poised enhancer. 
		os.system('bedtools intersect -a /mnt/silencer2/home/huh025/AL_HH/distal_h3k4me1_2_3_peak.bed -b /mnt/silencer2/home/huh025/AL_HH/method1_peak/h3k27ac.method1.peak -u > /mnt/silencer2/home/huh025/AL_HH/active_enhancer.bed')
		os.system('bedtools intersect -a /mnt/silencer2/home/huh025/AL_HH/distal_h3k4me1_2_3_peak.bed -b /mnt/silencer2/home/huh025/AL_HH/method1_peak/h3k27ac.method1.peak -v > /mnt/silencer2/home/huh025/AL_HH/poised_enhancer.bed')

		for i in all_input2:
			each=i.split()
			id=each[0]+'.'+each[1]+'.'+each[2]
			distal_dic[id]=1
		output2=open('/mnt/silencer2/home/huh025/AL_HH/active_enhancer.bed','w')
		output3=open('/mnt/silencer2/home/huh025/AL_HH/poised_enhancer.bed','w')
		print all_input1[0].split()[11]
		#find the index of h3k27ac
		for i in all_input1[1:]:
			each=i.split()
			id1=each[0]+'.'+each[1]+'.'+each[2]
			if float(each[11])>2:
				if distal_dic.has_key(id1):
					output2.write(each[0]+'\t'+each[1]+'\t'+each[2]+'\n')
			else:
				if distal_dic.has_key(id1):
					output3.write(each[0]+'\t'+each[1]+'\t'+each[2]+'\n')
 
			

def enhancer_associated_regions():

##check how the change of h3k4me1, h3k4me2, h3k4me3 in DKO cells associated with enhancer regions.

		input2=open('/mnt/silencer2/home/huh025/AL_HH/distal_h3k4me1_2_3_peak.bed','r')
		all_input2=input2.readlines()
		enhancer_dic={}
		for i in all_input2:
			each=i.split()
			id1=each[0]+'.'+each[1]+'.'+each[2]
			enhancer_dic[id1]=1
		input2.close()
		input3=open('/mnt/silencer2/home/huh025/AL_HH/all_peak_den_123_merged.txt','r')
		all_input3=input3.readlines()
		name=all_input3[0].split()
		i=0
		##find the index for h3k4me1,2,3
		while i<len(name):
			if name[i]=='h3k4me1':
				h3k4me1_idx=i
			if name[i]=='h3k4me1_dko':
				h3k4me1_dko_idx=i
			if name[i]=='h3k4me2':
				h3k4me2_idx=i
			if name[i]=='h3k4me2_dko':
                                h3k4me2_dko_idx=i
                        if name[i]=='h3k4me3':
                                h3k4me3_idx=i
                        if name[i]=='h3k4me3_dko':
                                h3k4me3_dko_idx=i
			i+=1
		each=all_input3[0].split()
		print each[h3k4me1_idx],each[h3k4me1_dko_idx],each[h3k4me2_idx],each[h3k4me2_dko_idx],each[h3k4me3_idx],each[h3k4me3_dko_idx]
		#find the toal decreased h3k4me3 regions, and the number of regions with concurrent decreased h3k4me1
		##decreased region are defined by over two fold enrichment compared to input, peak density less than half in DKO compared to WT.
		count1=0
		count2=0
		count3=0
		count4=0
		count5=0
		count6=0
		for i in all_input3[1:]:
			each=i.split()
			id1=each[0]+'.'+each[1]+'.'+each[2]
			if each[3].count('3')>1:
			  if float(each[h3k4me3_idx])>1:
				count6+=1
				if float(each[h3k4me3_dko_idx])/float(each[h3k4me3_idx])<0.5:
					if float(each[h3k4me1_idx])>1:
						if float(each[h3k4me1_dko_idx])/float(each[h3k4me1_idx])<0.5:
							#concurrent decreased regions
							count1+=1
							#count concurrent decreased regions by enhancer and others
							if enhancer_dic.has_key(id1):
								count4+=1
							else:
								count5+=1
						else:
							#other regions
							count2+=1
					else:
						count2+=1
					#total decreased h3k4me3 regions
					count3+=1
		print 'total h3k4me3'+'\t'+str(count6)+'\n'
		print 'total decreased h3k4me3'+'\t'+str(count3)+'\n'
		print 'concurrent decrease h3k4me3' +'\t'+ str(count1)+'\n'
		print 'concrrent decrease in enhancer'+'\t'+str(count4)+'\n'
		print 'concurrent decrease in others'+'\t'+str(count5)+'\n'
		print 'other regions'+'\t'+str(count2)+'\n'
		
		count1=0
        count2=0
        count3=0
        count4=0
        count5=0
        count6=0
        for i in all_input3[1:]:
                        each=i.split()
                        id1=each[0]+'.'+each[1]+'.'+each[2]
                        if each[3].count('2')>0:
                          if float(each[h3k4me2_idx])>1:
                                count6+=1
                                if float(each[h3k4me2_dko_idx])/float(each[h3k4me2_idx])<0.5:
                                        if float(each[h3k4me1_idx])>1:
                                                if float(each[h3k4me1_dko_idx])/float(each[h3k4me1_idx])<0.5:
                                                        #concurrent decreased regions
                                                        count1+=1
                                                        #count concurrent decreased regions by enhancer and others
                                                        if enhancer_dic.has_key(id1):
                                                                count4+=1
                                                        else:
                                                                count5+=1
                                                else:
                                                        #other regions
                                                        count2+=1
                                        else:
                                                count2+=1
                                        #total decreased h3k4me3 regions
                                        count3+=1
        print 'total h3k4me2'+'\t'+str(count6)+'\n'
        print 'total decreased h3k4me2'+'\t'+str(count3)+'\n'
        print 'concurrent decrease h3k4me2' +'\t'+ str(count1)+'\n'
		print 'concrrent decrease in enhancer'+'\t'+str(count4)+'\n'
        print 'concurrent decrease in others'+'\t'+str(count5)+'\n'
		print 'other regions'+'\t'+str(count2)+'\n'
                                                                              
		count1=0
		count2=0
		count3=0
        count4=0
        count5=0
        count6=0
		count7=0
		count8=0
		for i in all_input3[1:]:
                      	each=i.split()
                        id1=each[0]+'.'+each[1]+'.'+each[2]
			if each[3].count('1')>0:
                          if float(each[h3k4me1_idx])>1:
                                count6+=1
                                if float(each[h3k4me1_dko_idx])/float(each[h3k4me1_idx])<0.5:
                                        
                                                        #h3k4me1 decreased regions
                                                        count1+=1
                                                        #count decreased regions by enhancer and others
                                                        if enhancer_dic.has_key(id1):
                                                                count4+=1
                                                        else:
                                                                count5+=1
                                else:
                                                        #non-decreased regions
                                                        count2+=1
							#count decreased regions by enhancer and others
                                                        if enhancer_dic.has_key(id1):
                                                                count7+=1
                                                        else:
                                                                count8+=1
                                        
                                                
                                        #total decreased h3k4me3 regions
        print 'total h3k4me1'+'\t'+str(count6)+'\n'
        print 'total decreased h3k4me1'+'\t'+str(count1)+'\n'
        #print 'concurrent decrease h3k4me3' +'\t'+ str(count1)+'\n'
        print 'decreased in enhancer'+'\t'+str(count4)+'\n'
        print 'decreased in others'+'\t'+str(count5)+'\n'
        print 'other regions'+'\t'+str(count2)+'\n'
		print 'other region overlapped enhancer\t'+str(count7)+'\n'
        print 'other region overlapped non-enhancer\t'+str(count8)+'\n'                       
		count1=0
        count2=0
        count3=0
        count4=0
        count5=0
        count6=0
		count7=0
		count8=0
		for i in all_input3[1:]:
                        each=i.split()
                        id1=each[0]+'.'+each[1]+'.'+each[2]
			if each[3].count('2')>0:
                          if float(each[h3k4me2_idx])>1:
                                count6+=1
                                if float(each[h3k4me2_dko_idx])/float(each[h3k4me2_idx])<0.5:

                                                        #h3k4me2 decreased regions
                                                        count1+=1
                                                        #count decreased regions by enhancer and others
							if (each[h3k4me1_idx])>1:
                                                        	if enhancer_dic.has_key(id1):
                                                                	count4+=1
                                                        	else:
                                                                	count5+=1
							else:
								count5+=1
                                else:
                                                        #non-decreased regions
                                                        count2+=1
							#count decreased regions by enhancer and others
                                                        if (each[h3k4me1_idx])>1:
                                                                if enhancer_dic.has_key(id1):
                                                                        count7+=1
                                                                else:
                                                                        count8+=1
                                                        else:
                                                                count8+=1
		print 'total h3k4me2'+'\t'+str(count6)+'\n'
        print 'total decreased h3k4me2'+'\t'+str(count1)+'\n'
        #print 'concurrent decrease h3k4me3' +'\t'+ str(count1)+'\n'
        print 'decreased in enhancer'+'\t'+str(count4)+'\n'
        print 'decreased in others'+'\t'+str(count5)+'\n'
        print 'other regions'+'\t'+str(count2)+'\n'
		print 'other region overlapped enhancer\t'+str(count7)+'\n'
        print 'other region overlapped non-enhancer\t'+str(count8)+'\n'
		count1=0
        count2=0
        count3=0
        count4=0
        count5=0
        count6=0
		count7=0
		count8=0
        for i in all_input3[1:]:
                        each=i.split()
                        id1=each[0]+'.'+each[1]+'.'+each[2]
                        if each[3].count('3')>1:
                          if float(each[h3k4me3_idx])>1:
                                count6+=1
                                if float(each[h3k4me3_dko_idx])/float(each[h3k4me3_idx])<0.5:

                                                        #h3k4me2 decreased regions
                                                        count1+=1
                                                        #count decreased regions by enhancer and others
                                                        if (each[h3k4me1_idx])>1:
                                                                if enhancer_dic.has_key(id1):
                                                                        count4+=1
                                                                else:
                                                                        count5+=1
                                                        else:
                                                                count5+=1
                                else:
                                                        #non-decreased regions
                                                        count2+=1
							#count decreased regions by enhancer and others
                                                        if (each[h3k4me1_idx])>1:
                                                                if enhancer_dic.has_key(id1):
                                                                        count7+=1
                                                                else:
                                                                        count8+=1
                                                        else:
                                                                count8+=1
        print 'total h3k4me3'+'\t'+str(count6)+'\n'
        print 'total decreased h3k4me3'+'\t'+str(count1)+'\n'
		print 'decreased in enhancer'+'\t'+str(count4)+'\n'
        print 'decreased in others'+'\t'+str(count5)+'\n'
		print 'other regions'+'\t'+str(count2)+'\n'
		print 'other region overlapped enhancer\t'+str(count7)+'\n'
		print 'other region overlapped non-enhancer\t'+str(count8)+'\n'

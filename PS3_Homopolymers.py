# -*- coding: utf-8 -*-
"""
PS3 RNA-seq Data Exploration

Created on Wed Aug 24 20:14:40 2016

@author: jsmith


Find out approximately what proportion of the R1 reads contains a piece of a poly-A tail. 
Use regular expressions (with grep, awk, etc.) or a python script to accomplish this. 
Assume that anything composed of at least 15 consecutive “A”s is a poly-A tail.

Alsorunyourscripttocountoccurrencesof15ormoreconsecutive“C”s,“G”s,and “T”s and print 
their frequencies. 
"""

#Import Regular Expressions Module 
import re


#read in the fastq file 
file = "/home12/jsmith16/bi623/160825_PS3/PE_RNAseq_R1.fastq"
fastq = open(file, "r")

#compile the regular expression for each homopolymer nucleotide sequence 
#use {15,} to identify sequences with 15 or more of the nucleotides in a row 
polyA = re.compile("[A]{15,}")
polyT = re.compile("[T]{15,}")
polyC = re.compile("[C]{15,}")
polyG = re.compile("[G]{15,}")

#initialize variables to hold the count of each sequence read found with homopolymers
countA = 0
countT = 0 
countC = 0
countG = 0

#define a function to search each sequence read(the string) for poly nucleotides
def poly(string):
	global countA #make the variables in the function global, so that counts can be added 
	global countT #to the variables initialized outside the function
	global countC
	global countG
	if polyA.search(string): #search the read for 15 or more As
		countA += 1 #add one to the count for each read with a poly A sequence 
        #return countA #RETURN ends the function so dont use this here
	if polyT.search(string): #search for 15 or more Ts
		countT += 1  #add one to the count for each read with a poly T sequence 
	if polyC.search(string): #search for 15 or more Cs
		countC += 1  #add one to the count for each read with a poly C sequence 
	if polyG.search(string): #search for 15 or more Gs
		countG += 1   #add one to the count for each read with a poly G sequence  
    
#for loop to pick out the sequence reads with poly nucleotides    
i = 0 
for line in fastq:
    line = line.strip()
    if i % 4 == 1: #operation only on the second line of fastq, the sequence read
        poly(line) #call the function poly to search for the compiled regular expression   
    i += 1
    
#print the number of homopolymer sequences found in the fastq file
print("The number of polyA tails is", countA)
print("The number of polyC is", countC)
print("The number of polyT is", countT)
print("The number of polyG is", countG)


 
    
  
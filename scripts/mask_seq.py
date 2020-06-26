# this code was kindly provided by Nicola De Maio on 26/6/20

import sys
import os
import os.path
from os import path

alleles=["A","C","G","T"]
allelesLow=["a","c","g","t"]
allelesLowAmb=["a","c","g","t","n","-"]


infile = sys.argv[1]
outfile = sys.argv[2]

#read input alignment and create list of entry sequences.
def readFastaAlignmentKeepName(file,makeLower=False):
	samples=[]
	line=file.readline()
	while line!="" and line!="\n":
		if line[0]==">":
			name=line.replace("\n","").replace(">","")
			seq=""
			line=file.readline()
			while line!="" and line!="\n" and line[0]!=">":
				seq+=line.replace("\n","")
				line=file.readline()
			seq=seq.lower()
			samples.append([name,seq])
		else:
			print("problem with fasta format: line not recognised")
			print(line)
			exit()
	return samples
	

#alignment is assumed to be with respect to the reference (insertions with respect to the reference should have been removed)
infileh=open(infile)
samples=readFastaAlignmentKeepName(infileh)
infileh.close()
newSamples=[]
for s in samples:
	seq=s[1]
	newSeq=""
	
	#mask first 30 informative characters of the sequence
	p=0
	while seq[p]=="n" or seq[p]=="-":
		newSeq+=seq[p]
		p+=1
	for i in range(30):
		if seq[p]=="-":
			newSeq+="-"
		else:
			newSeq+="n"
		p+=1
		
	#mask last 30 positions of the sequence
	newSeqEnd=""
	p2=len(seq)-1
	while seq[p2]=="n" or seq[p2]=="-":
		newSeqEnd=seq[p2]+newSeqEnd
		p2-=1
	for i in range(30):
		if seq[p2]=="-":
			newSeqEnd="-"+newSeqEnd
		else:
			newSeqEnd="n"+newSeqEnd
		p2-=1
			
	#mask sites that have at least two ambiguities ("n" or "-") in a radius of 7bp around them.
	numN=0
	for i in range(15):
		if seq[p-7+i]=="n" or seq[p-7+i]=="-":
			numN+=1
	while p<=p2:
		if numN>1:
			if seq[p]=="-":
				newSeq+="-"
			else:
				newSeq+="n"
		else:
			newSeq+=seq[p]
		if (seq[p-7]=="n" or seq[p-7]=="-"):
			numN-=1
		p+=1
		if (seq[p+7]=="n" or seq[p+7]=="-"):
			numN+=1
	newSeq+=newSeqEnd
	if len(seq)!=len(newSeq):
		print('length error!')
		exit()
	newSamples.append([s[0],newSeq])
outfileh=open(outfile, 'w')
for s in newSamples:
	outfileh.write(">"+s[0]+"\n"+s[1]+"\n")
outfileh.close()

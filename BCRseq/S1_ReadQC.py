#!/usr/bin/env python
import os
import sys
import re
import time
import subprocess
import fnmatch


#######################################################################################
###  locate sample name and read/bam file
########################################################################################
currentPath=os.getcwd()
fastqFileFolder="/projects/b1042/Yuanqing/RNAseq/fastqFolder"
destF="/projects/b1042/Yuanqing/RNAseq/FQ_savefolder"
tempF='/projects/b1042/Yuanqing/RNAseq/BCR'
openFile=open("./pBCR_TLSpaper.txt","r").readlines()

#######################################################################################
### tools
########################################################################################
TrimTool="/projects/p31513/work/tools/Trimmomatic-0.39/trimmomatic-0.39.jar"


sampleName=[]
for line0 in openFile:
	line=line0.strip().split("\t")
	if line[0].startswith("ID"):
		continue
	else:
		nameH=line[0].replace("\x00","").replace("\r","")+"_HeavyChain"
		sampleName.append(nameH)
		nameL=line[0].replace("\x00","").replace("\r","")+"_LightChain"
		sampleName.append(nameL)
print (sampleName)


fastqcFolder=tempF+"/fastqc"
if not os.path.exists(fastqcFolder):
	os.makedirs(fastqcFolder)

fastqcFolderF=tempF+"/fastqcF"
if not os.path.exists(fastqcFolderF):
	os.makedirs(fastqcFolderF)

fastqFolder=destF+"/fastq"
if not os.path.exists(fastqFolder):
	os.makedirs(fastqFolder)

def fastqc ():
	for i in range(len(sampleName)):
		pbsfile=(tempF+"/rQC_%s.slurm" %(sampleName[i]))
		outfile=open(pbsfile,'w')
		outfile.write('#!/bin/bash\n')
		outfile.write('#SBATCH -A b1042\n')
		outfile.write('#SBATCH -p genomics\n')
		outfile.write('#SBATCH -N 1\n')
		outfile.write('#SBATCH --ntasks-per-node=2\n')
		outfile.write('#SBATCH -t 0:30:00\n')
		outfile.write('#SBATCH --mem-per-cpu=16G\n')
		outfile.write('#SBATCH -J QC_%s\n' %(sampleName[i]))
		outfile.write('#SBATCH -o %s/QC_%s.o%%j\n' %(tempF,sampleName[i]))
		outfile.write('#SBATCH -e %s/QC_%s.e%%j\n' %(tempF,sampleName[i]))

		outfile.write('#SBATCH --mail-user=yuanqing.yan@northwestern.edu\n')
		outfile.write('#SBATCH --mail-type=fail\n')


		outfile.write('date\n')
		outfile.write('cd %s\n' %(destF))
		outfile.write('module load fastqc/0.11.5\n')
		cp1=f"cp {fastqFileFolder}/{sampleName[i]}_R1.fastq.gz {tempF}/{sampleName[i]}_R1.fastq.gz\n"
		cp2=f"cp {fastqFileFolder}/{sampleName[i]}_R2.fastq.gz {tempF}/{sampleName[i]}_R2.fastq.gz\n"
		qcCom=f'fastqc {tempF}/{sampleName[i]}_R1.fastq.gz {tempF}/{sampleName[i]}_R2.fastq.gz -o {tempF}/fastqc\n'

		outfile.write("%s%s%s" %(cp1,cp2,qcCom))
		outfile.close()

		msubCmd=["sbatch %s" %pbsfile]
		subprocess.call(msubCmd,shell=True)


def trim():
	for i in range(len(sampleName)):
		pbsfile=(tempF+"/rT%s.slurm" %(sampleName[i]))
		outfile=open(pbsfile,'w')
		outfile.write('#!/bin/bash\n')
		outfile.write('#SBATCH -A b1042\n')
		outfile.write('#SBATCH -p genomics\n')
		outfile.write('#SBATCH -N 1\n')
		outfile.write('#SBATCH --ntasks-per-node=2\n')
		outfile.write('#SBATCH -t 0:30:00\n')
		outfile.write('#SBATCH --mem-per-cpu=16G\n')
		outfile.write('#SBATCH -J Trim_%s\n' %(sampleName[i]))
		outfile.write('#SBATCH -o %s/Trim_%s.o%%j\n' %(tempF,sampleName[i]))
		outfile.write('#SBATCH -e %s/Trim_%s.e%%j\n' %(tempF,sampleName[i]))

		outfile.write('#SBATCH --mail-user=yuanqing.yan@northwestern.edu\n')
		outfile.write('#SBATCH --mail-type=fail\n')


		outfile.write('date\n')
		outfile.write('cd %s\n' %(destF))
		outfile.write('module load fastqc/0.11.5\n')


		trim=f'java -jar {TrimTool} \
		PE \
		-phred33 \
		-threads 8 \
		{tempF}/{sampleName[i]}_R1.fastq.gz {tempF}/{sampleName[i]}_R2.fastq.gz \
		{destF}/fastq/{sampleName[i]}_R1_paired.fq.gz {destF}/fastq/{sampleName[i]}_R1_unpaired.fq.gz \
		{destF}/fastq/{sampleName[i]}_R2_paired.fq.gz {destF}/fastq/{sampleName[i]}_R2_unpaired.fq.gz \
		LEADING:20 \
		TRAILING:20 \
		SLIDINGWINDOW:4:20 \
		MINLEN:50\n'

		gzip1=f'gzip {destF}/fastq/{sampleName[i]}_R1_paired.fq.gz\n'
		gzip2=f'gzip {destF}/fastq/{sampleName[i]}_R2_paired.fq.gz\n'
		qcCom2=f'fastqc {destF}/fastq/{sampleName[i]}_R1_paired.fq.gz {destF}/fastq/{sampleName[i]}_R2_paired.fq.gz -o {tempF}/fastqcF\n'

		outfile.write("%s%s%s%s" %(trim,gzip1,gzip1,qcCom2))
		outfile.close()

		msubCmd=["sbatch %s" %pbsfile]
		subprocess.call(msubCmd,shell=True)



#fastqc ()
#trim()

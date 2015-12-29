__author__ = 'jb393'
import sys
import os
import re

vcf_file=open(sys.argv[1])
vcfTab=sys.argv[1].split("/")[-1]
snpVCF=open("SNP_"+vcfTab,"w")
indelVCF=open("INDEL_"+vcfTab,"w")

for data in vcf_file:
	start=re.compile("^\#")
	test2=start.search(data)
	if not test2:
 		data=data.strip()
		dataTab=data.split("\t")
		findSNP=re.compile("TSA=SNV")
		findDel=re.compile("TSA=deletion")
		findIns=re.compile("TSA=insertion")
		test=findSNP.search(dataTab[7])
		if test:
			snpVCF.write(data+"\n")
		else:
			indelVCF.write(data+"\n")

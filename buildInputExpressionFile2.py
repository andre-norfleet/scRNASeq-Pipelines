#!/usr/bin/env python
"""
author: jun ding
date: 2020-07-08

function: this code is used to generate  a tab-formated expression file required by scdiff software suite 

"""
import sys,os,pdb
import argparse 
import scipy.io
import csv
from scipy.sparse import csr_matrix
import gzip

def main():
	parser=argparse.ArgumentParser(description="convert 10x genomics mtx.matrix, barcodes.csv, and genes.tsv to a tab expression file that scdiff2 accepts")
	parser._action_groups.pop()
	required = parser.add_argument_group('required arguments')
	optional = parser.add_argument_group('optional arguments')
	required.add_argument('-i','--inputFolder',required=True,help='the input folder that stores mtx.matrix, barcodes.csv, and genes.tsv' +
																'you can also provide an optional meta.tsv to provide the time and label information for each of cells' +
																'if not provided, by default, all the cells will be set as time point 1 and label NA')
	args = parser.parse_args()
	exFn=args.inputFolder.strip()
	exFn=exFn.split("/")
	exFn=[item for item in exFn if item!=""]
	exFn="/".join(exFn)
	print(exFn)
	try:
		mat=scipy.io.mmread(f"{exFn}/matrix.mtx.gz")
		mat=mat.transpose()
		feature_ids = [row[0] for row in csv.reader(gzip.open(f"{exFn}/features.tsv.gz","rt"), delimiter="\t")]
		gene_names = [row[1] for row in csv.reader(gzip.open(f"{exFn}/features.tsv.gz","rt"), delimiter="\t")]
		barcodes = [row[0] for row in csv.reader(gzip.open(f"{exFn}/barcodes.tsv.gz","rt"), delimiter="\t")]
		mat3=scipy.io.mmread(f"C:/Users/Andre/Downloads/CPTfiltered_feature_bc_matrix/matrix.mtx.gz")
		mat3=mat3.transpose()
		feature_ids3 = [row[0] for row in csv.reader(gzip.open(f"C:/Users/Andre/Downloads/CPTfiltered_feature_bc_matrix/features.tsv.gz","rt"), delimiter="\t")]
		gene_names3 = [row[1] for row in csv.reader(gzip.open(f"C:/Users/Andre/Downloads/CPTfiltered_feature_bc_matrix/features.tsv.gz","rt"), delimiter="\t")]
		barcodes3 = [row[0] for row in csv.reader(gzip.open(f"C:/Users/Andre/Downloads/CPTfiltered_feature_bc_matrix/barcodes.tsv.gz","rt"), delimiter="\t")]
		mat1=scipy.io.mmread(f"C:/Users/Andre/Downloads/PPUfiltered_feature_bc_matrix/matrix.mtx.gz")
		mat1=mat1.transpose()
		feature_ids1 = [row[0] for row in csv.reader(gzip.open(f"C:/Users/Andre/Downloads/PPUfiltered_feature_bc_matrix/features.tsv.gz","rt"), delimiter="\t")]
		gene_names1 = [row[1] for row in csv.reader(gzip.open(f"C:/Users/Andre/Downloads/PPUfiltered_feature_bc_matrix/features.tsv.gz","rt"), delimiter="\t")]
		barcodes1 = [row[0] for row in csv.reader(gzip.open(f"C:/Users/Andre/Downloads/PPUfiltered_feature_bc_matrix/barcodes.tsv.gz","rt"), delimiter="\t")]
	except:
		mat=scipy.io.mmread(f"{exFn}/matrix.mtx")
		mat=mat.transpose()
		gene_names = [row[1] for row in csv.reader(open(f"{exFn}/features.tsv","rt"), delimiter="\t")]
		barcodes = [row[0] for row in csv.reader(open(f"{exFn}/barcodes.tsv","rt"), delimiter="\t")]
	
	mat=csr_matrix(mat)
	mat1=csr_matrix(mat1)
	mat3=csr_matrix(mat3)
	if os.path.exists(f"{exFn}/meta.tsv.gz"):
		meta=[row[1] for row in csv.read(gzip.open(f"{exFn}/meta.tsv.gz","rt"),delimiter="\t")]
	elif os.path.exists(f"{exFn}/meta.tsv"):
		meta=[row for row in csv.read(open(f"{exFn}/meta.tsv","rt"),delimiter="\t")]
	else:
		meta=[]
	
	FR=['cell','time','label']+gene_names
	#pdb.set_trace()
	#g=open(exFn+".E","a")
	g=open("C:/Users/Andre/Downloads/PPUCMUCPTfiltered_feature_bc_matrix4.E","a")
	g.write("\t".join(FR)+"\n")
	for i in range(len(barcodes)):
		ci=barcodes[i]
		#ci3=barcodes3[i]
		ti=meta[i][1] if len(meta)>0 else 2
		#ti3=meta[i][1] if len(meta)>0 else 3
		li=meta[i][2] if len(meta)>0 else 'NA'
		#li3=meta[i][3] if len(meta)>0 else 'CPT'
		ei=mat[i].toarray()[0]
		#ei3=mat3[i].toarray()[0]
		oi=[ci,ti,li]+list(ei)
		#oi3=[ci3,ti3,li3]+list(ei3)
		oi=[str(item) for item in oi]
		g.write("\t".join(oi)+"\n")
		#glook = "\t".join(oi)+"\n"
		print(i)
		#print(glook)
	#g.close()

	#FR3=['cell','time','label']+gene_names3
	#g3=open("C:/Users/Andre/Downloads/PPUCMUCPTfiltered_feature_bc_matrix.E","a")
	#g3.write("\t".join(FR3)+"\n")
	for i in range(len(barcodes3)):
		ci3=barcodes3[i]
		ti3=meta[i][1] if len(meta)>0 else 3
		li3=meta[i][2] if len(meta)>0 else 'NA'
		ei3=mat3[i].toarray()[0]
		oi3=[ci3,ti3,li3]+list(ei3)
		oi3=[str(item) for item in oi3]
		g.write("\t".join(oi3)+"\n")
		print(i)

	for i in range(len(barcodes1)):
		ci1=barcodes1[i]
		ti1=meta[i][1] if len(meta)>0 else 1
		li1=meta[i][2] if len(meta)>0 else 'NA'
		ei1=mat1[i].toarray()[0]
		oi1=[ci1,ti1,li1]+list(ei1)
		oi1=[str(item) for item in oi1]
		g.write("\t".join(oi1)+"\n")
		print(i)
	g.close()
	
if __name__=='__main__':
	main()
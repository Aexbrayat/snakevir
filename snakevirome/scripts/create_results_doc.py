#!/usr/local/bioinfo/python/2.7.9_build2/bin/python
# -*- coding: utf-8 -*-
import sys, os
from collections import Counter
import pandas as pd
import csv
from os.path import exists
import glob
# import pandas as pd
# import locale
# locale.setlocale(locale.LC_ALL, '')
path_files=os.listdir(sys.argv[1])
ext=sys.argv[2]
stats_id=sys.argv[3]
stats_seq=sys.argv[4]
out=sys.argv[5]

def WriteDictToCSV(csv_file,csv_columns,dict_data):
    with open(csv_file, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
        writer.writeheader()
        for data in dict_data:
            writer.writerow(data)
csv_columns=['Stats','All-sample']
dict_data = [{'Stats':"Total read R1"},{'Stats':"Total read R2"},{'Stats':"Avg read length R1"},{'Stats':"Avg read length R2"},{'Stats':"Reads 1 cleaned"},{'Stats':"Reads 2 cleaned"},
{'Stats':"Avg cleaned read length R1 (after cutadapt)"},{'Stats':"Avg cleaned read length R2 (after cutadapt)"},{'Stats':"Avg insert_size"},{'Stats':"Map on Diptera"},{'Stats':"UnMapped on Diptera"},{'Stats':"Map on bacteria"},{'Stats':"UnMapped on bacteria"},
{'Stats':"Total pairs"},{'Stats':"Combined pairs"},{'Stats':"Widows"},{'Stats':"Assembled"},{'Stats':"Reads spread across n contig (nb for each samples and Avg for all_samples)"},{'Stats':"Singletons"},{'Stats':"Nb total of contigs"},
{'Stats':"Min contigs length"},{'Stats':"Max contigs length"},{'Stats':"Avg contigs length"},{'Stats':"Contigs with viral hit"},{'Stats':"Nb of reads with viral hit 'inside' a contig"},{'Stats':"Nb of eads with viral hit as singleton"},{'Stats':"Reads with viral hit"},
{'Stats':"Nb of viral hit 10 =< Reads < 100 "},{'Stats':"Nb of viral hit 100 =< Reads < 1000"},{'Stats':"Nb of viral hit 1000 =< Reads < 10000"},{'Stats':"Nb of viral hit 10000=< Reads"},{'Stats':"Nb viral family"},{'Stats':"Nb viral genus"},{'Stats':"Nb viral species"}]
for i in range(27):
	dict_data[i]['All-sample']=int(0)
list_files=[]
for input_files in path_files :
	if input_files.endswith(ext):
		file_name=input_files.split(ext)[0]
		csv_columns.append(file_name)
		list_files.append(file_name)

for files_name in list_files :
	with open("logs/logscutadapt/"+files_name+"_1_rm_seq_adpt.log",'r') as cut1_1:
		for line in cut1_1:
			if "Total reads processed" in line:
				c = line.replace(',','').replace(' ','').replace('\n','').split(":")
				dict_data[0][files_name]=int(c[-1])
			if "Total written" in line:
				c1 = line.replace(',','').replace('b',':').replace('\n','').replace(' ','').split(":")
				dict_data[2][files_name]="%.2f" % (int(c1[1])/int(c[-1]))
	with open("logs/logscutadapt/"+files_name+"_2_rm_seq_adpt.log",'r') as cut1_2:
		for line in cut1_2:
			if "Total reads processed" in line:
				c = line.replace(',','').replace(' ','').replace('\n','').split(":")
				dict_data[1][files_name]=int(c[-1])
			if "Total written" in line:
				c1 = line.replace(',','').replace('b',':').replace('\n','').replace(' ','').split(":")
				dict_data[3][files_name]="%.2f" % (int(c1[1])/int(c[-1]))
	with open("logs/logscutadapt/"+files_name+"_1_qtrim.log",'r') as cut2_1:
		for line in cut2_1:
			if "Reads written" in line:
				c = line.replace(',','').replace('(',':').replace('\n','').replace(' ','').split(":")
				dict_data[4][files_name]=int(c[2])
			if "Total written" in line:
				c1 = line.replace(',','').replace('b',':').replace('\n','').replace(' ','').split(":")
				dict_data[6][files_name]="%.2f" % (int(c1[1])/int(c[2]))
	with open("logs/logscutadapt/"+files_name+"_2_qtrim.log",'r') as cut2_2:
		for line in cut2_2:
			if "Reads written" in line:
				c = line.replace(',','').replace('(',':').replace('\n','').replace(' ','').split(":")
				dict_data[5][files_name]=int(c[2])
			if "Total written" in line:
				c1 = line.replace(',','').replace('b',':').replace('\n','').replace(' ','').split(":")
				dict_data[7][files_name]="%.2f" % (int(c1[1])/int(c[2]))
	path="logs/insert_size/"+files_name+"_insert_size_metrics_*"
	files_insert=glob.glob(path)
	with open(files_insert[0],'r') as insert:
		for line in insert:
			if "LIBRARY	READ_GROUP" in line:
				line = next(insert)
				c = line.split("\t")
				dict_data[8][files_name]="%.2f" % (float(c[5]))

	with open("logs/logsFLASH/"+files_name+"_flash.log",'r') as flash:
		for line in flash:
			if "Total pairs" in line:
				c = line.replace('\n','').replace(' ','').split(":")
				dict_data[13][files_name]=int(c[1])
			if "Combined pairs" in line:
				c = line.replace('\n','').replace(' ','').split(":")
				dict_data[14][files_name]=int(c[1])
	path="logs/logs_coverage/"+files_name+"_coverage_*"
	files_cov=glob.glob(path)
	with open(files_cov[0],'r') as cov:
		for line in cov:
			if "Mapped" in line:
				c = line.split(":")
				dict_data[16][files_name]=int(c[1])
				line = next(cov)
				c = line.split(" ")
				dict_data[17][files_name]=int(c[0])
				line = next(cov)
				c = line.split(":")
				dict_data[18][files_name]=int(c[1])
	with open("logs/logs_contaminent/Stats_contaminent_"+files_name+".txt",'r') as conta:
		for line in conta:
			if "filter1_pair_R1" in line:
				c = line.replace('\n','').split(":")
				R1_h=c[1]
				line = next(conta)
				c = line.replace('\n','').split(":")
				R2_h=c[1]
				line = next(conta)
				c = line.replace('\n','').split(":")
				wi_h=c[1]
				line = next(conta)
				c = line.replace('\n','').split(":")
				R1_b=c[1]
				line = next(conta)
				c = line.replace('\n','').split(":")
				R2_b=c[1]
				line = next(conta)
				c = line.replace('\n','').split(":")
				wi_b=c[1]
				dict_data[15][files_name]=int(wi_b)
				dict_data[10][files_name]=int(R1_h)+int(R2_h)+int(wi_h)
				dict_data[12][files_name]=int(R1_b)+int(R2_b)+int(wi_b)
				dict_data[9][files_name]=(int(dict_data[4][files_name])+int(dict_data[5][files_name]))-int(dict_data[10][files_name])
				dict_data[11][files_name]=dict_data[10][files_name]-dict_data[12][files_name]
		for i in (0,1,4,5,9,10,11,12,13,14,15,16,18):
			dict_data[i]['All-sample']+=int(dict_data[i][files_name])
		for i in (2,3,6,7,8,17):
			if dict_data[i]['All-sample']==0:
				dict_data[i]['All-sample']=round((float(dict_data[i][files_name])), 2)
			else:
				dict_data[i]['All-sample']= round((float(dict_data[i]['All-sample']+float(dict_data[i][files_name]))/2), 2)
#
	df = pd.read_csv(stats_id, delimiter='\t')
	stat_vir=df[files_name] != 0
	df1=(df[stat_vir].family.value_counts())
	dict_data[31][files_name]=df1.count()
	df1_1=(df.family.value_counts())
	dict_data[31]['All-sample']=df1_1.count()
	df2=(df[stat_vir].genus.value_counts())
	dict_data[32][files_name]=df2.count()
	df2_2=(df.genus.value_counts())
	dict_data[32]['All-sample']=df2_2.count()
	df3=(df[stat_vir].species.value_counts())
	dict_data[33][files_name]=df3.count()
	df3_3=(df.species.value_counts())
	dict_data[33]['All-sample']=df3_3.count()
	raws = df[[files_name]]
	sum0=sum10=sum100=sum1000=sum10000=sum0c=0
	numcont=0
	for index, raw in raws.iterrows():
		numcont+=1
		if raw[files_name] != 0 :
			sum0+=raw[files_name]
		if 10 <= raw[files_name] < 100 :
			sum10+=1
		if 100 <= raw[files_name] <1000 :
			sum100+=1
		if 1000 <= raw[files_name] < 10000 :
			sum1000+=1
		if raw[files_name] >= 10000 :
			sum10000+=1
	dict_data[26][files_name]=sum0
	dict_data[27][files_name]=sum10
	dict_data[28][files_name]=sum100
	dict_data[29][files_name]=sum1000
	dict_data[30][files_name]=sum10000
	df = pd.read_csv(stats_seq, delimiter='\t')
	raws = df[[files_name,'qseqid']]
	in_contigs=0
	numcont=0
	for index, raw in raws.iterrows():
		numcont+=1
		if raw[files_name] != 0 :
			if raw['qseqid'].startswith('k') or raw['qseqid'].startswith('C'):
				in_contigs+=raw[files_name]
	out_oncitgs=sum0-in_contigs
	dict_data[23]['All-sample']=numcont
	dict_data[24][files_name]=in_contigs
	dict_data[25][files_name]=out_oncitgs
	for i in range(19,24):
		dict_data[i][files_name]='-'
	for i in range(24,27):
		dict_data[i]['All-sample']+=dict_data[i][files_name]


tot_hit10=tot_hit100=tot_hit1000=tot_hit10000=0
df = pd.read_csv(stats_id, delimiter='\t')
for index, raw in df.iterrows():
	sum_raw=0
	for input_files in path_files :
		if input_files.endswith(ext):
			files_name=input_files.split(ext)[0]
			sum_raw+=raw[files_name]
	if 10 <= sum_raw< 100 :
		tot_hit10+=1
	if 100 <= sum_raw < 1000 :
		tot_hit100+=1
	if 1000 <= sum_raw < 10000 :
		tot_hit1000+=1
	if sum_raw >= 10000 :
		tot_hit10000+=1

dict_data[27]['All-sample']=tot_hit10
dict_data[28]['All-sample']=tot_hit100
dict_data[29]['All-sample']=tot_hit1000
dict_data[30]['All-sample']=tot_hit10000

path=os.listdir("logs/logsAssembly/")
for input_files in path :
	if input_files.endswith("Assembly.log"):
		with open("logs/logsAssembly/"+input_files,'r') as asb:
			data=pd.read_csv(asb, sep="\t", header=None)
			dict_data[19]['All-sample']=len(data.index)
			c=data.min()
			dict_data[20]['All-sample']=c.values[1]
			c=data.max()
			dict_data[21]['All-sample']=c.values[1]
			c=data.mean()
			dict_data[22]['All-sample']="%.2f" %(c.values[0])

csv_file = out
WriteDictToCSV(csv_file,csv_columns,dict_data)

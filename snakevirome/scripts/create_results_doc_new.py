import sys, os
from collections import Counter
import pandas as pd
import csv
from os.path import exists
import glob
import numpy as np
# import pandas as pd
# import locale
# locale.setlocale(locale.LC_ALL, '')
path_files=os.listdir(sys.argv[1])
ext1=sys.argv[2]
ext2=sys.argv[3]
ext_file=sys.argv[4]
run=sys.argv[5]
out=sys.argv[9]


stats_df=pd.read_csv(sys.argv[6], sep=',').drop_duplicates()
lin_df=pd.read_csv(sys.argv[7], sep=',').drop_duplicates()
count_df_r=pd.read_csv(sys.argv[8], sep=',').drop_duplicates()

stats_df=stats_df[(stats_df.tax_id.isin(lin_df.tax_id))]
count_df=count_df_r[(count_df_r.qseqid.isin(stats_df.qseqid))]
stats_df=stats_df[(stats_df.qseqid.isin(count_df.qseqid))]



df_s_l=stats_df.merge(lin_df, on='tax_id', how='left')
df=df_s_l.merge(count_df, on='qseqid', how='left')
df = df[df.iloc[: , -1].notna()]
tmp_df=df.iloc[:, np.r_[0,30:len(df. columns)]]
tmp_df.to_csv(r'tmp_df.csv')


df.to_csv(r'df.csv')
csv_columns=['All-sample']
list_files=[]
for input_files in path_files :
	if input_files.endswith(ext1+ext_file):
		file_name=input_files.split(ext1+ext_file)[0]
		csv_columns.append(file_name)
		list_files.append(file_name)

result_table = pd.DataFrame(columns=csv_columns,index=["Total read R1","Total read R2","Avg read length R1","Avg read length R2","Reads 1 cleaned","Reads 2 cleaned","Avg cleaned read length R1 (after cutadapt)","Avg cleaned read length R2 (after cutadapt)","Avg insert_size","Map on Diptera","UnMapped on Diptera","Map on bacteria","UnMapped on bacteria","Total pairs","Combined pairs","Widows","Assembled","Duplicates","Reads_count","Singletons","Nb of contigs","Min contigs length","Max contigs length","Avg contigs length","Contigs with viral hit","Nb of reads with viral hit 'inside' a contig","Nb of reads with viral hit as singleton","Reads with viral hit","% viral reads","Nb of viral hit 1 =< Reads < 10 ","Nb of viral hit 10 =< Reads < 100 ","Nb of viral hit 100 =< Reads < 1000","Nb of viral hit 1000 =< Reads < 10000","Nb of viral hit 10000=< Reads","Nb viral family","Nb viral genus","Nb viral species"])

contig_length=pd.read_csv("logs/logsAssembly/"+run+"_assembly_stats.txt", sep='	')
print(contig_length)
stats_contigs=contig_length.iloc[:, 1].describe()

contig_length.columns.values[0] = "qseqid"


n=0
for files_name in list_files :
	with open("logs/logscutadapt/"+files_name+ext1+"_cut1.log",'r') as cut1_1:
		for line in cut1_1:
			if "Total reads processed" in line:
				c = line.replace(',','').replace(' ','').replace('\n','').split(":")
				result_table.at['Total read R1', files_name] =int(c[-1])
			if "Total written" in line:
				c1 = line.replace(',','').replace('b',':').replace('\n','').replace(' ','').split(":")
				result_table.at['Avg read length R1', files_name]="%.2f" % (int(c1[1])/int(c[-1]))
	cut1_1.close()
	with open("logs/logscutadapt/"+files_name+ext2+"_cut1.log",'r') as cut1_2:
		for line in cut1_2:
			if "Total reads processed" in line:
				c = line.replace(',','').replace(' ','').replace('\n','').split(":")
				result_table.at['Total read R2', files_name]=int(c[-1])
			if "Total written" in line:
				c1 = line.replace(',','').replace('b',':').replace('\n','').replace(' ','').split(":")
				result_table.at['Avg read length R2', files_name]="%.2f" % (int(c1[1])/int(c[-1]))
	cut1_2.close()
	with open("logs/logscutadapt/"+files_name+ext1+"_cut2.log",'r') as cut2_1:
		for line in cut2_1:
			if "Reads written" in line:
				c = line.replace(',','').replace('(',':').replace('\n','').replace(' ','').split(":")
				result_table.at['Reads 1 cleaned', files_name]=int(c[2])
			if "Total written" in line:
				c1 = line.replace(',','').replace('b',':').replace('\n','').replace(' ','').split(":")
				result_table.at['Avg cleaned read length R1 (after cutadapt)', files_name]="%.2f" % (int(c1[1])/int(c[2]))
	cut2_1.close()
	with open("logs/logscutadapt/"+files_name+ext2+"_cut2.log",'r') as cut2_2:
		for line in cut2_2:
			if "Reads written" in line:
				c = line.replace(',','').replace('(',':').replace('\n','').replace(' ','').split(":")
				result_table.at['Reads 2 cleaned', files_name]=int(c[2])
			if "Total written" in line:
				c1 = line.replace(',','').replace('b',':').replace('\n','').replace(' ','').split(":")
				result_table.at['Avg cleaned read length R2 (after cutadapt)', files_name]="%.2f" % (int(c1[1])/int(c[2]))
	cut2_2.close()
	path="logs/insert_size/"+files_name+"_insert_size_metrics_"+run+".txt"
	files_insert=glob.glob(path)
	with open(files_insert[0],'r') as insert:
		for line in insert:
			if "LIBRARY	READ_GROUP" in line:
				line = next(insert)
				c = line.split("\t")
				result_table.at['Avg insert_size', files_name]="%.2f" % (float(c[5]))
	insert.close()
	with open("logs/logsFLASH/"+files_name+"_flash.log",'r') as flash:
		for line in flash:
			if "Total pairs" in line:
				c = line.replace('\n','').replace(' ','').split(":")
				result_table.at['Total pairs', files_name]=int(c[1])
			if "Combined pairs" in line:
				c = line.replace('\n','').replace(' ','').split(":")
				result_table.at['Combined pairs', files_name]=int(c[1])
	flash.close()
	
	path="logs/logs_coverage/"+files_name+"_coverage*"
	files_cov=glob.glob(path)
	with open(files_cov[0],'r') as cov:
		for line in cov:
			if "Mapped" in line:
				c = line.split(":")
				result_table.at['Reads_count', files_name]=int(c[1])
				print(files_name,result_table.at['Reads_count', files_name])
				line = next(cov)
				c = line.split(" ")
				result_table.at['Nb of contigs', files_name]=int(c[0])
				line = next(cov)
				c = line.split(":")
				result_table.at['Singletons', files_name]=int(c[1])
	cov.close()
		
	with open("logs/logsDuplicates/"+files_name+"_duplicates_pairs_"+run+".txt",'r') as dp,  open("logs/logsDuplicates/"+files_name+"_duplicates_wi_"+run+".txt",'r') as dw :
		for line in dp :
			if "READ:" in line:
				reads_pair = int(line.replace('\n','').split(": ")[-1])
			if "DUPLICATE TOTAL:" in line:
				dup_pair = int(line.replace('\n','').split(": ")[-1])
				
		for line in dw:
			if "READ:" in line:
				reads_wi = int(line.replace('\n','').split(": ")[-1])
			if "DUPLICATE TOTAL:" in line:
				dup_wi = int(line.replace('\n','').split(": ")[-1])
		result_table.at['Duplicates', files_name]=dup_pair+dup_wi
		result_table.at['Assembled', files_name]=result_table.at['Duplicates', files_name]+result_table.at['Reads_count', files_name]
	

	
	with open("logs/logs_contaminent/Stats_contaminent_"+files_name+".txt",'r') as conta:

		for line in conta:
			if "Host_pair_R1" in line:
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
				result_table.at['Widows', files_name]=int(wi_b)
				result_table.at['UnMapped on Diptera', files_name]=int(R1_h)+int(R2_h)+int(wi_h)
				result_table.at['UnMapped on bacteria', files_name]=int(R1_b)+int(R2_b)+int(wi_b)
				result_table.at['Map on Diptera', files_name]=(int(result_table.at['Reads 1 cleaned', files_name])+int(result_table.at['Reads 2 cleaned', files_name]))-int(result_table.at['UnMapped on Diptera', files_name])
				result_table.at['Map on bacteria', files_name]=result_table.at['UnMapped on Diptera', files_name]-result_table.at['UnMapped on bacteria', files_name]
	conta.close()
	
	result_table.at['Nb viral family', files_name]=df[["family", files_name]][df[files_name]> 1].family.value_counts().count()
	result_table.at['Nb viral genus', files_name]=df[["genus", files_name]][df[files_name]> 1].genus.value_counts().count()
	result_table.at['Nb viral species', files_name]=df[["species", files_name]][df[files_name]> 1].species.value_counts().count()
	raws = df[["species",files_name]].groupby('species',as_index=False).agg({files_name: 'sum' }).drop(columns='species')
	
	ranges = [1,10,100,1000,10000]
	result_table.at['Reads with viral hit', files_name]=raws.values.sum()
	result_table.at['Nb of viral hit 1 =< Reads < 10 ', files_name]=raws.groupby(pd.cut(raws[files_name], ranges)).count().iloc[0][files_name]
	result_table.at['Nb of viral hit 10 =< Reads < 100 ', files_name]=raws.groupby(pd.cut(raws[files_name], ranges)).count().iloc[1][files_name]
	result_table.at['Nb of viral hit 100 =< Reads < 1000', files_name]=raws.groupby(pd.cut(raws[files_name], ranges)).count().iloc[2][files_name]
	result_table.at['Nb of viral hit 1000 =< Reads < 10000', files_name]=raws.groupby(pd.cut(raws[files_name], ranges)).count().iloc[3][files_name]
	result_table.at['Nb of viral hit 10000=< Reads', files_name]=raws[raws[files_name]> 10000].count().values[0]
	result_table.at['Contigs with viral hit', files_name]=count_df[count_df[files_name]> 0].count().values[0]
	result_table.at['Min contigs length', files_name]=int(contig_length[(contig_length.qseqid.isin(count_df_r[count_df_r[files_name]!= 0][["qseqid",files_name]] .qseqid))].iloc[:, 1].describe().iloc[3])

	result_table.at['Max contigs length', files_name]=int(contig_length[(contig_length.qseqid.isin(count_df_r[count_df_r[files_name]!= 0][["qseqid",files_name]] .qseqid))].iloc[:, 1].describe().iloc[7])

	result_table.at['Avg contigs length', files_name]=round(contig_length[(contig_length.qseqid.isin(count_df_r[count_df_r[files_name]!= 0][["qseqid",files_name]] .qseqid))].iloc[:, 1].describe().iloc[1],2)
	result_table.at["Nb of reads with viral hit 'inside' a contig", files_name]=count_df[count_df[["qseqid",files_name]].qseqid.astype(str).str.startswith('k') | count_df[["qseqid",files_name]].qseqid.astype(str).str.startswith('C')][files_name].sum()
	result_table.at["Nb of reads with viral hit as singleton", files_name]=count_df[~count_df[["qseqid",files_name]].qseqid.astype(str).str.startswith('k') & ~count_df[["qseqid",files_name]].qseqid.astype(str).str.startswith('C')][files_name].sum()
	result_table.at['% viral reads', files_name]=round(result_table.loc["Reads with viral hit", files_name] *100 / (result_table.loc["Total read R1", files_name]+result_table.loc["Total read R1", files_name]),2)
result_table.at['Nb of contigs', 'All-sample']=int(count_df_r.shape[0])-1
result_table.at['Min contigs length', 'All-sample']=int(stats_contigs.iloc[3])
result_table.at['Max contigs length', 'All-sample']=int(stats_contigs.iloc[7])
result_table.at['Avg contigs length', 'All-sample']=round(stats_contigs.iloc[1],2)
result_table.at['Contigs with viral hit', 'All-sample']=count_df.shape[0]
result_table.at['Reads with viral hit', 'All-sample']=count_df.iloc[:, 1:].values.sum()
result_table.at['Nb viral family', 'All-sample']=df.family.value_counts().count()
result_table.at['Nb viral genus', 'All-sample']=df.genus.value_counts().count()
result_table.at['Nb viral species', 'All-sample']=df.species.value_counts().count()


spp=df.iloc[:,np.r_[df.columns.get_loc("species"), df_s_l.shape[1]:df.shape[1]]].groupby('species',as_index=False).agg('sum').drop(columns='species').sum(axis=1)
result_table.at['Nb of viral hit 1 =< Reads < 10 ',  'All-sample']=spp.groupby(pd.cut(spp, ranges)).count().iloc[0]
result_table.at['Nb of viral hit 10 =< Reads < 100 ',  'All-sample']=spp.groupby(pd.cut(spp, ranges)).count().iloc[1]
result_table.at['Nb of viral hit 100 =< Reads < 1000',  'All-sample']=spp.groupby(pd.cut(spp, ranges)).count().iloc[2]
result_table.at['Nb of viral hit 1000 =< Reads < 10000',  'All-sample']=spp.groupby(pd.cut(spp, ranges)).count().iloc[3]
result_table.at['Nb of viral hit 10000=< Reads',  'All-sample']=spp[spp > 10000].count()

row_sum=["Total read R1","Total read R2","Reads 1 cleaned","Reads 2 cleaned","Map on Diptera","UnMapped on Diptera","Map on bacteria","UnMapped on bacteria","Total pairs","Combined pairs","Widows","Assembled","Duplicates","Reads_count","Singletons","Nb of reads with viral hit 'inside' a contig","Nb of reads with viral hit as singleton","Reads with viral hit"]

row_mean=["Avg read length R1","Avg read length R2","Avg cleaned read length R1 (after cutadapt)","Avg cleaned read length R2 (after cutadapt)","Avg insert_size"]
print("ok3")

for row in row_sum:
	result_table.at[row,'All-sample']=result_table.loc[row,result_table.columns != 'All-sample'].apply(int).sum(axis=0)
for row in row_mean:
	result_table.at[row,'All-sample']=result_table.loc[row,result_table.columns != 'All-sample'].apply(float).mean(axis=0).round(2)


result_table.at['% viral reads',  'All-sample']=round(result_table.loc["Reads with viral hit", 'All-sample'] *100 / (result_table.loc["Total read R1", 'All-sample']+result_table.loc["Total read R1", 'All-sample']),2)


result_table.to_csv(r'result_table_dup.csv')



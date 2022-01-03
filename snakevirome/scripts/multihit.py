import sys, os
from ete3 import NCBITaxa
from collections import Counter
import pandas as pd
import csv
import subprocess
import numpy as np
from collections import Counter

blast_contig_hit_file=sys.argv[1]
blast_contig_tax_file=sys.argv[2]
tax_table_file=sys.argv[3]
out_file=sys.argv[4]
#Taxonomy_full_file=sys.argv[2]
#Ref_taxonomy_file=sys.argv[3]
#ICTV_taxo_file=sys.argv[4]
#Host_taxo_file=sys.argv[5]
#output_stat= sys.argv[6]
#output_seq= sys.argv[7]

if not os.path.exists('tmp'):
	os.makedirs('tmp')

blast_contig_hit=pd.read_csv(blast_contig_hit_file, sep=',')
blast_contig_hit=blast_contig_hit.fillna("none")
print("1")

blast_contig_tax=pd.read_csv(blast_contig_tax_file, sep=',')
blast_contig_tax=blast_contig_tax.fillna("none")

tax_table=pd.read_csv(tax_table_file, sep=',')
tax_table=tax_table[tax_table.superkingdom == "Viruses"]
tax_table=tax_table.fillna("NA")


start_table=blast_contig_hit.merge(blast_contig_tax, on="tax_id", how='left')

start_table=start_table.fillna("none")
out_table=start_table.drop_duplicates(subset=['qseqid'])
out_table=out_table[out_table.species != "none"]
start_table=start_table[start_table['qseqid'].isin(out_table['qseqid'])]
out_table=out_table.drop_duplicates(subset=['tax_id'])



out_table.to_csv(r'out_table.csv', sep=',')
start_table.to_csv(r'start_table.csv', sep=',')

def Average(lst):
    return sum(lst) / len(lst)
def GetKey(dictA,val):
   for key, value in dictA.items():
      if val == value:
         return key
      continue
print("2")
qseq_dic={}
tax_dic={}
pident={}
n=0
for qseqid_0 in  start_table["qseqid"].unique():
	n+=1

	if qseqid_0 not in qseq_dic:
		if start_table[start_table ["qseqid"]== qseqid_0]["species"].tolist()[0]!="none":
			qseq_dic[qseqid_0]=start_table[start_table ["qseqid"]== qseqid_0]["family"].tolist()
			qseq_dic[qseqid_0].insert(0,start_table[start_table ["qseqid"]== qseqid_0]["tax_id"].tolist()[0])
print("ok1")
for qseqid_0 in  start_table["qseqid"].unique():
	taxid_0=start_table[start_table ["qseqid"]== qseqid_0]["tax_id"].tolist()[0]
	if taxid_0 not in pident:
		pident[taxid_0]={}
	for index, row in start_table[start_table ["qseqid"]== qseqid_0].iterrows():
		if row["family"] not in pident[taxid_0]:
			pident[taxid_0][row["family"]]=[row["pident"]]
		else:
			pident[taxid_0][row["family"]].append(row["pident"])

for k, v in pident.items():
	for k1 in pident[k]:
		if "none" not in pident[k][k1]:
			pident[k][k1]=round(Average(pident[k][k1]),2)
print("ok2")
df=pd.DataFrame.from_dict(qseq_dic,orient='index')
df=df.rename(columns={0: "tax_id"})
df.fillna("none")
df1 = df.groupby("tax_id").agg(lambda x: x.tolist())
df1.fillna("none")

df1['combined']=df1[1]+df1[2]+df1[3]+df1[4]+df1[5]
print("ok3")
df1.fillna("none")
df1.to_csv(r'df1.csv', sep=',')
for index, row in df1.iterrows():
	row.combined=["none" if v is None else v for v in row.combined]
	if index not in tax_dic:
		tax_dic[index]=[len(row.combined)/5]
	if len(set(row.combined))==1:
		tax_dic[index].extend([len(set(row.combined)),dict(Counter(row.combined))])
	else:
		row.combined[:]=(value for value in row.combined if value != "none")
		tax_dic[index].extend([len(set(row.combined)),dict(Counter(row.combined))])
	if tax_dic[index][2]:
		tax_dic[index].extend([sum(list(tax_dic[index][2].values())),max(list(tax_dic[index][2].values()))])
		tax_dic[index].append(pident[index])
		try :
			tax_dic[index].append(max({k: tax_dic[index][5][k] for k in tax_dic[index][5].keys() if k not in ('none')}, key={k: tax_dic[index][5][k] for k in tax_dic[index][5].keys() if k not in ('none')}.get))
		except ValueError:
			tax_dic[index].append(max( tax_dic[index][5], key=tax_dic[index][5].get))
		try :
			tax_dic[index].append(max(list({k: tax_dic[index][5][k] for k in tax_dic[index][5].keys() if k not in ('none')}.values())))
		except ValueError:
			tax_dic[index].append(max(list(tax_dic[index][5].values())))
		tax_dic[index].append(max(tax_dic[index][2].values())*100/sum(list(tax_dic[index][2].values())))
		if max(tax_dic[index][2], key=tax_dic[index][2].get)==tax_dic[index][6]:
			if max(list(tax_dic[index][2].values()))*100/sum(list(tax_dic[index][2].values())) >= 75 :
					tax_dic[index].append(max(tax_dic[index][2], key=tax_dic[index][2].get))
			elif tax_dic[index][1]==3 and max(list(tax_dic[index][2].values()))*100/sum(list(tax_dic[index][2].values())) >= 60:
				tax_dic[index].append(max(tax_dic[index][2], key=tax_dic[index][2].get))
			elif tax_dic[index][1]==2 and max(list(tax_dic[index][2].values()))*100/sum(list(tax_dic[index][2].values())) > 50:
				tax_dic[index].append(max(tax_dic[index][2], key=tax_dic[index][2].get))

print("ok4")
for k, v in tax_dic.items():
	if v[1] == 0:
		del tax_dic[k]



df2=pd.DataFrame.from_dict(tax_dic,orient='index')
df2new = df2[[9]]
df2new=df2new.fillna("NA")
df2new.index.name = 'tax_id'

#tax_table.set_index('tax_id').join(df2new.set_index('tax_id'))
df2new.reset_index(inplace=True)

tax_table['tax_id'] = tax_table['tax_id'].astype(int)

tax_table=tax_table.merge(df2new, on="tax_id", how='left')

mask = (tax_table['family'] == "NA") & (tax_table[9] != "NA") & (tax_table[9] != "none")
tax_table['family'][mask] = tax_table[9][mask]
df2.columns=['n_contigs','n_species','fam_out','n_out','max_out','fam_pid','best_predict_fam','max_pid','per','final_fam']

df3 = df2[['fam_out','fam_pid','final_fam']]
print(df3)
df3.index.name = 'tax_id'
df3.reset_index(inplace=True)
print("ok5")

tax_table=tax_table.merge(df3, on="tax_id", how='left')

#data['column2'][mask] = 3

#if tax_table['family']=="" :
##	print(tax_table['family'])
#def fx(x):
#    if tax_table['family']=="NA" and tax_table[9]!="none" and tax_table[9]!="NA":
#        return tax_table[9]
#    else:
#        return tax_table['family']

#tax_table['family']=tax_table.apply(lambda x : fx(x),axis=1)
del tax_table[9]

#tax_table.to_csv(r'tax_table.csv', sep='\t', index=False)
tax_table.to_csv(sys.argv[4], sep='\t', index=False, mode = 'w')
df3.to_csv("Taxonomy/comp_inf.csv", sep='\t', index=False, mode = 'w')
#tax_table.to_csv(r out_file, sep=',')
##df.to_csv(r'qseq_dic.csv', sep=',')
##df1.to_csv(r'qseq_dic_taxid.csv', sep=',')
#df2.to_csv(r'df2.csv', sep=',')

##		tax_dic[index].extend([sum(tax_dic[index][2].values()),max(tax_dic[index][2].values()),max(tax_dic[index][2], key=tax_dic[index][2].get)])
#	max_key = max(tax_dic[index][2], key=tax_dic[index][2].get)

#	print(index, len(set(row.combined)),dict(Counter(row.combined)))

#df2 = df.groupby("tax_id").agg(lambda x: ','.join(x).split(','))
#df1=df1.merge(df2, on="tax_id", how='left')
#df['combined']= df.to_str(values).tolist()
#print(df1.groupby(0))

#df2.to_csv(r'tax_dic.csv', sep=',')
#qseq_dic.to_csv("qseq_dic.csv")
#		qseq_dic[qseqid_0]=[blast_contig_hit.loc[blast_contig_hit["qseqid"]== qseqid_0]["qseqid"].count(),blast_contig_hit[blast_contig_hit ["qseqid"]== qseqid_0]["family"].tolist()]
#	n+=1
#	print(n,blast_contig_hit[blast_contig_hit["qseqid"]== qseqid]["qseqid"])
#	print(qseqid,blast_contig_hit.loc[blast_contig_hit["qseqid"]== qseqid].counts(),blast_contig_hit[blast_contig_hit ["qseqid"]== qseqid]["tax_id"].nunique())
#	print(qseqid,blast_contig_hit.loc[blast_contig_hit["qseqid"]== qseqid]["qseqid"].count(),blast_contig_hit[blast_contig_hit ["qseqid"]== qseqid]["family"].nunique())
#	df1=blast_contig_hit[blast_contig_hit ["qseqid"]== qseqid]["tax_id"].value_counts()
#for hits in blast_contig_hit:
#	
#	ref_lineage.append(lineage.replace('\n','').split(','))
#print(blast_contig_hit)


#ref_lineage=[]
#ICTV_lineage=[]
#host_lineage=[]
#host_dict={}
#Taxonomy=open(Taxonomy_file, 'r')
#Taxonomy_full=open(Taxonomy_full_file, 'r')
#Ref_taxonomy=open(Ref_taxonomy_file, 'r')
#ICTV_taxo=open(ICTV_taxo_file, 'r')
#host_taxo=open(Host_taxo_file, 'r')

#for lineage in Ref_taxonomy:
#	ref_lineage.append(lineage.replace('\n','').split(','))


#for lineage in host_taxo:
#	host_lineage=(lineage.replace('\n','').split(','))
#	if host_lineage[0] != "virus_tax_id":
#		if host_lineage[0] not in host_dict:
#			host_dict[host_lineage[0]]=[host_lineage[1:]]
#		else:
#			host_dict[host_lineage[0]].append(host_lineage[1:])
#for lineage in ICTV_taxo:
#	ICTV_lineage.append(lineage.replace('\n','').split(','))

#liste_full=[]
#liste_fullseqid=[]

#for lineage in Taxonomy_full:
#	rank=lineage.replace('\n','').replace('\r','').split(",")
#	if rank[0]=="qseqid":
#		ranks_full=rank[0:16]+["clade"]+rank[16:21]+["cluster"]+["Nucleic_acid"]+rank[21:22]+["hosts_gb"]+["hosts_tax_id"]+["hosts_superkingdom"]+["hosts_class"]+["hosts_order"]+["hosts_family"]+["hosts_subfamily"]+["hosts_genus"]+["hosts_species"]+rank[22:]
#	elif rank[14]=="Viruses":
#		ranks_full=rank[0:16]+["NA"]+rank[16:21]+["NA"]+["NA"]+rank[21:22]+["NA"]+["NA"]+["NA"]+["NA"]+["NA"]+["NA"]+["NA"]+["NA"]+["NA"]+rank[22:]
#	else:
#		continue
#	liste_full.append(ranks_full)

#host_all=[]
#for ranks_full in liste_full:
#	index=liste_full.index(ranks_full)
#	if liste_full[index][14]=="Viruses":
#		for ictv_rank in ICTV_lineage:
#			if liste_full[index][21]==ictv_rank[3] and  liste_full[index][21]!="NA" :
#				if liste_full[index][20]=="NA":
#					liste_full[index][20]=ictv_rank[2]
#				if liste_full[index][18]=="NA":
#					liste_full[index][18]=ictv_rank[1]
#				if liste_full[index][17]=="NA":
#					liste_full[index][17]=ictv_rank[0]
#				if liste_full[index][23]=="NA":
#					liste_full[index][23]=ictv_rank[4]

#		for ref_rank1 in ref_lineage:

#			if liste_full[index][21]==ref_rank1[6] and  liste_full[index][21]!="NA" :
#				if liste_full[index][20]=="NA":
#					liste_full[index][20]=ref_rank1[3]
#				if liste_full[index][18]=="NA":
#					liste_full[index][18]=ref_rank1[2]
#				if liste_full[index][17]=="NA":
#					liste_full[index][17]=ref_rank1[1]
#				if liste_full[index][23]=="NA":
#					liste_full[index][23]=ref_rank1[5]
#				if liste_full[index][22]=="NA":
#					liste_full[index][22]=ref_rank1[4]
#				if liste_full[index][16]=="NA":
#					liste_full[index][16]=ref_rank1[0]

#		if liste_full[index][18]!="NA" and liste_full[index][22]=="NA":
#			liste_full[index][22]=liste_full[index][18]
#		if liste_full[index][19]!="NA" and liste_full[index][22]=="NA":
#			liste_full[index][22]=liste_full[index][16]
#		if liste_full[index][16]!="NA" and liste_full[index][22]=="NA":
#			liste_full[index][22]=liste_full[index][17]

#		if liste_full[index][23]=="NA":
#			if "ssRNA negative-strand viruses" in liste_full[index][24]:
#				liste_full[index][23]="ssRNA(-)"
#			if "ssRNA positive-strand viruses" in liste_full[index][24]:
#				liste_full[index][23]="ssRNA(+)"
#			if "dsRNA viruses" in liste_full[index][24]:
#				liste_full[index][23]="dsRNA"
#			if "RNA satellites" in liste_full[index][24]:
#				liste_full[index][23]="RNA satellites"
#			if "ssDNA viruses" in liste_full[index][24]:
#				liste_full[index][23]="ssDNA(+/-)"
#			if "dsDNA viruses" in liste_full[index][24]:
#				liste_full[index][23]="dsDNA"
#			if "unclassified ssRNA viruses" in liste_full[index][24]:
#				liste_full[index][23]="ssRNA"
#			if "unclassified RNA viruses" in liste_full[index][24]:
#				 liste_full[index][23]="RNA"
#		if  liste_full[index][1] not in seqid :
#			seqid.append(liste_full[index][1])
#		if ranks_full[13] in host_dict:
#			print(ranks_full[13])
#			host_id=""
#			superkingdom=[]
#			clas=[]
#			order=[]
#			family=[]
#			subfamily=[]
#			genus=[]
#			species=[]
#			for hosts in host_dict[ranks_full[13] ]:
#				if host_id=="" and hosts[0]!='':
#					host_id+=str(hosts[0]).replace("[","").replace("]","")
#				else:
#					host_id+="|"+str(hosts[0])
#				if hosts[1] not in superkingdom:
#					superkingdom.append(hosts[1])
#				if hosts[2] not in clas:
#					clas.append(hosts[2])
#				if hosts[3] not in order:
#					order.append(hosts[3])
#				if hosts[4] not in family:
#					family.append(hosts[4])
#				if hosts[5] not in subfamily:
#					subfamily.append(hosts[5])
#				if hosts[6] not in genus:
#					genus.append(hosts[6])
#				if hosts[7] not in species:
#					species.append(hosts[7])
#			liste_full[index][26]=(str(host_id).replace("']","").replace("['",""))
#			liste_full[index][27]=(str(superkingdom).replace("]","").replace("[","").replace("'","").replace(", ","|"))
#			liste_full[index][28]=(str(clas).replace("]","").replace("[","").replace("'","").replace(", ","|"))
#			liste_full[index][29]=(str(order).replace("]","").replace("[","").replace("'","").replace(", ","|"))
#			liste_full[index][30]=(str(family).replace("]","").replace("[","").replace("'","").replace(", ","|"))
#			liste_full[index][31]=(str(subfamily).replace("]","").replace("[","").replace("'","").replace(", ","|"))
#			liste_full[index][32]=(str(genus).replace("]","").replace("[","").replace("'","").replace(", ","|"))
#			liste_full[index][33]=(str(species).replace("]","").replace("[","").replace("'","").replace(", ","|"))


#ncbi = NCBITaxa()
#liste_host=[]
#stat_by_seq=[]
#n=0

#print("ok2")
#desired_ranks = ['superkingdom', 'class', 'order', 'family', 'subfamily','genus','species']
#for ids in seqid:
#	print(n,"/",len(seqid))
#	n+=1
#	ranks2lineage={}
#	liste=[]
#	fieldnames=[]
#	fieldnames1=[]
#	host_list=[]
#	host_id="NA"
##	command = 'esearch -db protein -query '+ str(ids) +'| efetch -format gb | grep -oP "/host=\K.*"'
##	proc = subprocess.Popen(command,stdout=subprocess.PIPE,shell=True)
##	(out, err) = proc.communicate()
##	host_gb=out.decode("utf-8").rstrip().replace('"','').replace("[u'","")
##	cleaned_host_gb=host_gb.split(' (', 1)[0]
#	cleaned_host_gb=''
#	if cleaned_host_gb=="mosquito" or cleaned_host_gb=="mosquitoes":
#		cleaned_host_gb="Culicoidea"
#	if cleaned_host_gb!='':
#		liste.extend([ids , cleaned_host_gb])
#		liste_host.append(liste)
#		name2taxid = ncbi.get_name_translator([cleaned_host_gb])
#		fieldnames.append(host_gb)
#		if cleaned_host_gb in name2taxid:
#			host_id=int(name2taxid[cleaned_host_gb][0])
#			lineage = ncbi.get_lineage(host_id)
#			fieldnames.append(host_id)
#	if cleaned_host_gb!='' and cleaned_host_gb in name2taxid:
#		lineage2name = ncbi.translate_to_names(lineage)
#		lineage2ranks=ncbi.get_rank(lineage)
#		Dic_lineage2name=dict(zip(lineage, lineage2name))
#		for (taxid,rank) in lineage2ranks.items():
#			if rank in desired_ranks:
#				if not rank in ranks2lineage:
#					ranks2lineage[rank]=[taxid]
#				else:
#					ranks2lineage[rank].append(taxid)
#		ranks2names={}
#		for rank in ranks2lineage:
#			ranks2names[rank]=""
#			for i in range(len(ranks2lineage[rank])):
#				ranks2names[rank]+=Dic_lineage2name[ranks2lineage[rank][i]]+";"
#			ranks2names[rank]=ranks2names[rank][:-1]
#		ranks2names['host_id']=host_id
#		for rank in desired_ranks:
#			fieldnames1.append(ranks2names.get(rank, 'NA'))
#	host_list=[ids]+fieldnames+fieldnames1
#	host_all.append(host_list)


#print("ok3")
#dict_host={}
#for seq in liste_full:
#	if seq[0] == "qseqid":
#		stat_by_seq.append(seq[0:])
#	else:
#		for host in host_all:
#			l=[]
#			if seq[1]==host[0]:
#				if len(host)==2:
#					l=seq[0:25]+host[1:]+seq[26:]
#				elif len(host)>2:
#					l=seq[0:25]+host[1:]+seq[34:]
#				else:
#					l=seq[0:]
#				stat_by_seq.append(l)
#				if len(host)>=2:
#					for idx_host, val_host in enumerate(host):
#						if  isinstance(val_host,str):
#							val_host=val_host
#						if seq[13] not in dict_host :
#							dict_host[seq[13]]={idx_host:[val_host]}
#						elif idx_host not in dict_host[seq[13]]:
#							dict_host[seq[13]].update({idx_host:[val_host]})
#						elif  val_host!='NA' and val_host not in dict_host[seq[13]][idx_host]:
#							dict_host[seq[13]][idx_host].append(val_host)


#for ranks_full in liste_full:
#	index=liste_full.index(ranks_full)
#	if str(liste_full[index][13]) in dict_host and len(dict_host[str(liste_full[index][13])])>2:
#		liste_full[index][25]=str(dict_host[str(liste_full[index][13])][2]).replace("]","").replace("[","").replace("'","").replace(", ","|")
#		liste_full[index][26]=str(dict_host[str(liste_full[index][13])][3]).replace("]","").replace("[","").replace("'","").replace(", ","|")
#		liste_full[index][27]=str(dict_host[str(liste_full[index][13])][4]).replace("]","").replace("[","").replace("'","").replace(", ","|")
#		liste_full[index][28]=str(dict_host[str(liste_full[index][13])][5]).replace("]","").replace("[","").replace("'","").replace(", ","|")
#		liste_full[index][29]=str(dict_host[str(liste_full[index][13])][6]).replace("]","").replace("[","").replace("'","").replace(", ","|")
#		liste_full[index][30]=str(dict_host[str(liste_full[index][13])][7]).replace("]","").replace("[","").replace("'","").replace(", ","|")
#		liste_full[index][31]=str(dict_host[str(liste_full[index][13])][8]).replace("]","").replace("[","").replace("'","").replace(", ","|")
#		liste_full[index][32]=str(dict_host[str(liste_full[index][13])][9]).replace("]","").replace("[","").replace("'","").replace(", ","|")

#liste_species=[]
#for lineage in Taxonomy:
#	rank=lineage.replace('\n','').replace('\r','').split(",")
#	if rank[0]=="tax_id":
#		ranks=rank[0:3]+["clade"]+rank[3:8]+["cluster"]+["Nucleic_acid"]+["hosts_gb"]+["hosts_tax_id"]+["hosts_superkingdom"]+["hosts_class"]+["hosts_order"]+["hosts_family"]+["hosts_subfamily"]+["hosts_genus"]+["hosts_species"]+rank[9:]
#	elif rank[1]=="Viruses":
#		ranks=rank[0:3]+["NA"]+rank[3:8]+["NA"]+["NA"]+["NA"]+["NA"]+["NA"]+["NA"]+["NA"]+["NA"]+["NA"]+["NA"]+["NA"]+rank[9:]
#	else:
#		continue
#	liste_species.append(ranks)


#for ranks in liste_species:
#	index=liste_species.index(ranks)
#	if liste_species[index][1]=="Viruses":
#		for ictv_rank in ICTV_lineage:
#			if liste_species[index][8]==ictv_rank[3] and liste_species[index][8]!="NA":
#				if liste_species[index][7]=="NA":
#					liste_species[index][7]=ictv_rank[2]
#				if liste_species[index][5]=="NA":
#					liste_species[index][5]=ictv_rank[1]
#				if liste_species[index][4]=="NA":
#					liste_species[index][4]=ictv_rank[0]
#				if liste_species[index][10]=="NA":
#					liste_species[index][10]=ictv_rank[4]


#		for ref_rank1 in ref_lineage:
#			if liste_species[index][8]==ref_rank1[6] and liste_species[index][8]!="NA" :
#				if liste_species[index][7]=="NA":
#					liste_species[index][7]=ref_rank1[3]
#				if liste_species[index][5]=="NA":
#					liste_species[index][5]=ref_rank1[2]
#				if liste_species[index][4]=="NA":
#					liste_species[index][4]=ref_rank1[1]
#				if liste_species[index][10]=="NA":
#					liste_species[index][10]=ref_rank1[5]
#				if liste_species[index][9]=="NA":
#					liste_species[index][9]=ref_rank1[4]

#		if liste_species[index][5]!="NA" and liste_species[index][9]=="NA":
#			liste_species[index][9]=liste_species[index][5]
#		if liste_species[index][7]!="NA" and liste_species[index][9]=="NA":
#			liste_species[index][9]=liste_species[index][7]
#		if liste_species[index][4]!="NA" and liste_species[index][9]=="NA":
#			liste_species[index][9]=liste_species[index][4]

#		if liste_species[index][10]=="NA":
#			if "ssRNA negative-strand viruses" in liste_species[index][11]:
#				liste_species[index][10]="ssRNA(-)"
#			if "ssRNA positive-strand viruses" in liste_species[index][11]:
#				liste_species[index][10]="ssRNA(+)"
#			if "dsRNA viruses" in liste_species[index][11]:
#				liste_species[index][10]="dsRNA"
#			if "RNA satellites" in liste_species[index][11]:
#				liste_species[index][10]="RNA satellites"
#			if "ssDNA viruses" in liste_species[index][11]:
#				liste_species[index][10]="ssDNA(+/-)"
#			if "dsDNA viruses" in liste_species[index][11]:
#				liste_species[index][10]="dsDNA"
#			if "unclassified ssRNA viruses" in liste_species[index][11]:
#				liste_species[index][10]="ssRNA"
#			if "unclassified RNA viruses" in liste_species[index][11]:
#				liste_species[index][10]="RNA"
#		if ranks[0] in host_dict:
#			host_id=""
#			superkingdom=[]
#			clas=[]
#			order=[]
#			family=[]
#			subfamily=[]
#			genus=[]
#			species=[]
#			for hosts in host_dict[ranks[0]]:
#				if host_id=="" and hosts[0]!='':
#					host_id+=str(hosts[0]).replace("[","").replace("]","")
#				else:
#					host_id+="|"+str(hosts[0])
#				if hosts[1] not in superkingdom:
#					superkingdom.append(hosts[1])
#				if hosts[2] not in clas:
#					clas.append(hosts[2])
#				if hosts[3] not in order:
#					order.append(hosts[3])
#				if hosts[4] not in family:
#					family.append(hosts[4])
#				if hosts[5] not in subfamily:
#					subfamily.append(hosts[5])
#				if hosts[6] not in genus:
#					genus.append(hosts[6])
#				if hosts[7] not in species:
#					species.append(hosts[7])
#			liste_species[index][12]=(str(host_id).replace("']","").replace("['",""))
#			liste_species[index][13]=(str(superkingdom).replace("]","").replace("[","").replace("'","").replace(", ","|"))
#			liste_species[index][14]=(str(clas).replace("]","").replace("[","").replace("'","").replace(", ","|"))
#			liste_species[index][15]=(str(order).replace("]","").replace("[","").replace("'","").replace(", ","|"))
#			liste_species[index][16]=(str(family).replace("]","").replace("[","").replace("'","").replace(", ","|"))
#			liste_species[index][17]=(str(subfamily).replace("]","").replace("[","").replace("'","").replace(", ","|"))
#			liste_species[index][18]=(str(genus).replace("]","").replace("[","").replace("'","").replace(", ","|"))
#			liste_species[index][19]=(str(species).replace("]","").replace("[","").replace("'","").replace(", ","|"))
#			if str(liste_species[index][0]) in dict_host:
#				liste_species[index][11]=str(dict_host[str(liste_species[index][0])][1]).replace("]","").replace("[","").replace("'","").replace(", ","|")
#		else:
#			if str(liste_species[index][0]) in dict_host and len(dict_host[str(liste_species[index][0])])>2:
#				liste_species[index][12]=str(dict_host[str(liste_species[index][0])][2]).replace("]","").replace("[","").replace("'","").replace(", ","|")
#				liste_species[index][13]=str(dict_host[str(liste_species[index][0])][3]).replace("]","").replace("[","").replace("'","").replace(", ","|")
#				liste_species[index][14]=str(dict_host[str(liste_species[index][0])][4]).replace("]","").replace("[","").replace("'","").replace(", ","|")
#				liste_species[index][15]=str(dict_host[str(liste_species[index][0])][5]).replace("]","").replace("[","").replace("'","").replace(", ","|")
#				liste_species[index][16]=str(dict_host[str(liste_species[index][0])][6]).replace("]","").replace("[","").replace("'","").replace(", ","|")
#				liste_species[index][17]=str(dict_host[str(liste_species[index][0])][7]).replace("]","").replace("[","").replace("'","").replace(", ","|")
#				liste_species[index][18]=str(dict_host[str(liste_species[index][0])][8]).replace("]","").replace("[","").replace("'","").replace(", ","|")
#				liste_species[index][19]=str(dict_host[str(liste_species[index][0])][9]).replace("]","").replace("[","").replace("'","").replace(", ","|")
#			elif str(liste_species[index][0]) in dict_host and len(dict_host[str(liste_species[index][0])])>=2:
#				liste_species[index][11]=str(dict_host[str(liste_species[index][0])][1]).replace("]","").replace("[","").replace("'","").replace(", ","|")
#			else:
#				continue
#list_vir_ids=[]
#dic_host_ids={}
#to_write=[]

#for listes in liste_species:
#	if listes[0]=="tax_id":
#		to_write.append(listes[0:11]+["host_taxon"]+listes[11:])
#	else:
#		to_write.append(listes[0:11]+["NA"]+listes[11:])
#	if "|" in listes[12] :
#		dic_host_ids[listes[0]]=(listes[12].split('|'))
#	else :
#		dic_host_ids[listes[0]]=([listes[12]])

#ranks2lineage={}
#ncbi = NCBITaxa()
#dic={}
#test=[]
#for key in dic_host_ids:
#	for ids in dic_host_ids[key]:
#		if ids != "NA" and ids.isdigit():
#			try:
#				lineage=ncbi.get_lineage(ids)
#			except:
#				pass
#			lineage2name = ncbi.translate_to_names(lineage)
#			lineage2ranks=ncbi.get_rank(lineage)
#			Dic_lineage2name=dict(zip(lineage, lineage2name))
#			for (taxid,rank) in lineage2ranks.items():
#				if taxid in Dic_lineage2name:
#					if key not in dic:
#						dic[key]={}
#					if rank not in dic[key] :
#						dic[key][rank]=[Dic_lineage2name[taxid]]
#					elif Dic_lineage2name[taxid] not in dic[key][rank]:
#						dic[key][rank].append(Dic_lineage2name[taxid])

#	if key in dic:
#		if 'kingdom' in dic[key] and len(dic[key]['kingdom'])==1 and dic[key]['kingdom'][0]!='NA' and dic[key]['kingdom'][0]!='Metazoa':
#			for line in to_write :
#				if key in line:
#					line[11]=str(dic[key]['kingdom'][0])

#		if 'superkingdom' in dic[key] and len(dic[key]['superkingdom'])==1 and dic[key]['superkingdom'][0]!="Eukaryota":
#			for line in to_write:
#				if key in line:
#					line[11]=str(dic[key]['superkingdom'][0])
#		if  'order' in dic[key] and 'class' in dic[key] and len(dic[key]['order'])==1 and len(dic[key]['class'])==1:
#			if dic[key]['order'][0]=='Diptera':
#				if 'family' in dic[key] and dic[key]['family'][0]=='Culicidae':
#					for line in to_write:
#						if key in line:
#							line[11]=str(dic[key]['family'][0])
#				else:
#					for line in to_write:
#						if key in line:
#							line[11]=str(dic[key]['order'][0])
#		for line in to_write:
#			if key in line and  line[11]=='NA' and 'class' in dic[key] and len(dic[key]['class'])==1:
#				if dic[key]['class'][0]=='Insecta':
#							line[11]=str(dic[key]['phylum'][0])
#				elif dic[key]['class'][0]!='Mammalia':
#					if 'Vertebrata' in dic[key]['no rank']:
#								line[11]='Vertebrata'
#				elif dic[key]['class'][0]=='Mammalia':
#							line[11]=str(dic[key]['class'][0])
#			if key in line and  'phylum' in dic[key] and  line[11]=='NA' and dic[key]['phylum'][0]=='Arthropoda':
#						line[11]='Arthropoda'
#			if key in line and  line[11]=='NA':
#				if len(line[13].split('|'))==1:
#					line[11]='Other'


#with open(output_stat, "w") as f:
#    writer = csv.writer(f)
#    writer.writerows(to_write)
#with open(output_seq, "w") as g:
#    writer = csv.writer(g)
#    writer.writerows(stat_by_seq)
#f.close
#g.close

import sys, os
from ete3 import NCBITaxa
from collections import Counter
import pandas as pd
import csv
import subprocess
import pickle
from Bio import Entrez
Entrez.email = "exbrayat.antoni@gmail.com"
Entrez.api_key = "a92c4a6007614fcd3fe29e002a41cf08c008"
Taxonomy_file=sys.argv[1]
seq_hit_file=sys.argv[2]
Host_taxo_file=sys.argv[3]
output= sys.argv[4]

Taxonomy=open(Taxonomy_file, 'r')
seq_hit=open(seq_hit_file, 'r')
host_taxo=open(Host_taxo_file, 'r')

host_lineage=[]
host_dict={}


for lineage in host_taxo:
	host_lineage=(lineage.replace('\n','').split(','))

	if host_lineage[0] != "virus_tax_id":
		if host_lineage[0] not in host_dict:
			host_dict[host_lineage[0]]=[host_lineage[1:2]]
		else :
			host_dict[host_lineage[0]].append(host_lineage[1:2])

tax_ids=[]
otu_tax={}
taxo_dic=pd.read_csv(sys.argv[1], index_col=0).to_dict()


for lineage in Taxonomy:
	rank=lineage.replace('\n','').replace('\r','').split(",")
	if rank[0]!="tax_id":
		if rank[0] not in tax_ids:
			tax_ids.append(rank[0])
		if rank[0] not in otu_tax and  rank[0]!="tax_id": 
			otu_tax[rank[0]]={"hosts_taxon":[],"hosts_gb":[],"hosts_tax_id":[],"hosts_superkingdom":[],"hosts_kingdom":[],
						"hosts_phylum":[],"hosts_class":[],"hosts_order":[],"hosts_family":[],"hosts_subfamily":[],"hosts_genus":[],"hosts_species":[]}
						
seq_tax={}
for seq in seq_hit:
	seq=seq.replace('\n','').replace('\r','').split(",")
	if seq[-1] in tax_ids:
		if seq[1] not in seq_tax: 
			seq_tax[seq[1]]={"tax_id":seq[-1],"hosts_taxon":"","hosts_gb":[],"hosts_tax_id":[],"hosts_superkingdom":"","hosts_kingdom":"",
						"hosts_phylum":"","hosts_class":"","hosts_order":"","hosts_family":"","hosts_subfamily":"","hosts_genus":"","hosts_species":""}

ncbi = NCBITaxa()
liste_host=[]
out=[]
n=0
host_all=[]
dic={}


for sequence in seq_tax:
	seq=seq_tax[sequence]["tax_id"]
	if seq in host_dict:
		host_id=[]
		superkingdom=[]
		clas=[]
		order=[]
		family=[]
		subfamily=[]
		genus=[]
		species=[]
		for hosts in host_dict[seq]:
			if  hosts[0]!='' and hosts[0] not in host_id:
				host_id.append(hosts[0])
		if host_id !='':
			seq_tax[sequence]["hosts_tax_id"]=host_id
	print(n,"/",len(seq_tax))
	n+=1
	ranks2lineage={}
	liste=[]
	fieldnames=[]
	fieldnames1=[]
	host_list=[]
	host_id="NA"
	taxid="txid"+ str(seq)
	print(taxid)
	handle = Entrez.esearch(db="protein", retmax=20, term= str(taxid), idtype="acc")
	record = Entrez.read(handle)
	sseqids=record["IdList"]
	if sequence in sseqids:
		sseqids.remove(sequence.replace('"',''))
	sseqids.insert(0, sequence.replace('"',''))
	cleaned_host_gb=""
	source=""
	for ids in sseqids:
		handle = Entrez.efetch(db="protein", id=str(ids), retmode="xml")
		records = Entrez.read(handle)
		records=records[0]['GBSeq_feature-table'][0]['GBFeature_quals']
		for record in records:
			if "host" in record.values(): 
				cleaned_host_gb=record['GBQualifier_value']
				species =  cleaned_host_gb.replace(" ", "+").strip()
				search = Entrez.esearch(term=species, db="taxonomy", retmode="xml")
				org_txi = Entrez.read(search)
				print("host ", cleaned_host_gb)
				if not org_txi['IdList']:
					continue
				else:
					host_id = org_txi['IdList'][0]
			elif "host" not in record.values() and "isolation_source" in record.values():
				source=record['GBQualifier_value']
		if source !="":
			break		
		if cleaned_host_gb !="":
			break
	if cleaned_host_gb=='' and source !="": 
		cleaned_host_gb=source

	if cleaned_host_gb!='':
		seq_tax[sequence]["hosts_gb"].append(cleaned_host_gb)
		liste.extend([sequence , cleaned_host_gb])
		liste_host.append(liste)
		
	if host_id not in seq_tax[sequence]["hosts_tax_id"] and host_id!="NA":
		seq_tax[sequence]["hosts_tax_id"].append(host_id)
	for hosts_ids in seq_tax[sequence]["hosts_tax_id"]:
		if hosts_ids != "" and len(seq_tax[sequence]["hosts_tax_id"])>0:
			try:
				lineage=ncbi.get_lineage(hosts_ids)
			except:
				continue
			lineage2name = ncbi.translate_to_names(lineage)
			lineage2ranks=ncbi.get_rank(lineage)
			Dic_lineage2name=dict(zip(lineage, lineage2name))
			for (taxid,rank) in lineage2ranks.items():
				if taxid in Dic_lineage2name:
					if sequence not in dic:
						dic[sequence]={}
					if rank not in dic[sequence] :
						dic[sequence][rank]=[Dic_lineage2name[taxid]]
					elif Dic_lineage2name[taxid] not in dic[sequence][rank]:
						dic[sequence][rank].append(Dic_lineage2name[taxid])
			if sequence in dic:
				if dic[sequence].get('superkingdom')!=None:
					seq_tax[sequence]["hosts_superkingdom"]=(dic[sequence].get('superkingdom'))
				if dic[sequence].get('kingdom')!=None:
					seq_tax[sequence]["hosts_kingdom"]=(dic[sequence].get('kingdom'))
				if dic[sequence].get('phylum')!=None:
					seq_tax[sequence]["hosts_phylum"]=(dic[sequence].get('phylum'))
				if dic[sequence].get('class')!=None:
					seq_tax[sequence]["hosts_class"]=(dic[sequence].get('class'))
				if dic[sequence].get('order')!=None:
					seq_tax[sequence]["hosts_order"]=(dic[sequence].get('order'))
				if dic[sequence].get('family')!=None:
					seq_tax[sequence]["hosts_family"]=(dic[sequence].get('family'))
				if dic[sequence].get('subfamily')!=None:
					seq_tax[sequence]["hosts_subfamily"]=(dic[sequence].get('subfamily'))
				if dic[sequence].get('genus')!=None:
					seq_tax[sequence]["hosts_genus"]=(dic[sequence].get('genus'))
				if dic[sequence].get('species')!=None:
					seq_tax[sequence]["hosts_species"]=(dic[sequence].get('species'))
					
#a_file = open("host.pkl", "rb")
#seq_tax = pickle.load(a_file)



a_file = open("host2.pkl", "wb")
pickle.dump(seq_tax, a_file)
a_file.close()

###a_file = open("host.pkl", "rb")
###output = pickle.load(a_file)
###print(output)


for key in seq_tax:

	if "Bacteria" in str(seq_tax[key]["hosts_superkingdom"])and not seq_tax[key]["hosts_taxon"]:
		seq_tax[key]["hosts_taxon"]=["Bacteria"]
#	if "Arthropoda"  in seq_tax[key]["hosts_phylum"] and "Chordata"  in seq_tax[key]["hosts_phylum"]and not seq_tax[key]["hosts_taxon"]:
#		seq_tax[key]["hosts_taxon"]=["Arbovirus"]
	elif len(seq_tax[key]["hosts_class"])==1 :
		if "Mammalia" in str(seq_tax[key]["hosts_class"])and not seq_tax[key]["hosts_taxon"]:
			seq_tax[key]["hosts_taxon"]=["Vertebrate"]
		elif "Aves" in str(seq_tax[key]["hosts_class"])and not seq_tax[key]["hosts_taxon"]:
			seq_tax[key]["hosts_taxon"]=["Vertebrate"]
		elif "Arthropoda" in str(seq_tax[key]["hosts_phylum"]):

			if "Insecta" in str(seq_tax[key]["hosts_class"]):
				if "Diptera" in str(seq_tax[key]["hosts_order"]) and  not seq_tax[key]["hosts_taxon"]:
					seq_tax[key]["hosts_taxon"]=["Diptera"]
				else:
					seq_tax[key]["hosts_taxon"]=["Insecta_other"]
			else:
				seq_tax[key]["hosts_taxon"]=["Arthropoda_other"]
		elif "Viridiplantae" in str(seq_tax[key]["hosts_kingdom"]):
			seq_tax[key]["hosts_taxon"]=["Viridiplantae"]
		elif "Fungi" in str(seq_tax[key]["hosts_kingdom"]):
			seq_tax[key]["hosts_taxon"]=["Fungi"]

	elif len(seq_tax[key]["hosts_class"])>1 :
		if "Mammalia" in str(seq_tax[key]["hosts_class"])and not seq_tax[key]["hosts_taxon"]:
			seq_tax[key]["hosts_taxon"]=["Vertebrate"]
		elif "Aves" in str(seq_tax[key]["hosts_class"]) and  not seq_tax[key]["hosts_taxon"]:
			seq_tax[key]["hosts_taxon"]=["Vertebrate"]
		elif "Arthropoda" in str(seq_tax[key]["hosts_phylum"]):
			if "Insecta" in str(seq_tax[key]["hosts_class"]):
				if "Diptera" in str(seq_tax[key]["hosts_order"]) and  not not seq_tax[key]["hosts_taxon"]:
					seq_tax[key]["hosts_taxon"]=["Diptera"]
				else:
					seq_tax[key]["hosts_taxon"]=["Insecta_other"]
			elif not seq_tax[key]["hosts_taxon"]:
				seq_tax[key]["hosts_taxon"]=["Arthropoda_other"]
		elif "Viridiplantae" in str(seq_tax[key]["hosts_kingdom"])  and not seq_tax[key]["hosts_taxon"]:
			seq_tax[key]["hosts_taxon"]=["Viridiplantae"]
		elif  not seq_tax[key]["hosts_taxon"] and seq_tax[key]["hosts_tax_id"][0]!="NA" :
			seq_tax[key]["hosts_taxon"]=["Other"]
	elif len(seq_tax[key]["hosts_tax_id"])>=1 and not seq_tax[key]["hosts_taxon"] and seq_tax[key]["hosts_tax_id"][0]!="NA":
		seq_tax[key]["hosts_taxon"]=["Other"]

	elif not seq_tax[key]["hosts_taxon"]:
		if "phage" in taxo_dic["species"][int(seq_tax[key]["tax_id"])]:
			seq_tax[key]["hosts_taxon"]=["Bacteria"]


			
	for key1, value in seq_tax[key].items():
		tax_id_m=str(seq_tax[key]["tax_id"])
		if key1 !="tax_id":
			otu_tax[tax_id_m][key1]=list(set(otu_tax[tax_id_m][key1]) | set(seq_tax[key][key1]))
		if isinstance(value, list):
			seq_tax[key].update({key1:";".join(value)})

for key in otu_tax:
	for key1, value in otu_tax[key].items():
		if isinstance(value, list):
			otu_tax[key].update({key1:";".join(value)})


#df=pd.DataFrame.from_dict(seq_tax,orient='index')
#df.index.name = 'sseqid'
#df.reset_index(inplace=True)
#column_names=["sseqid","tax_id","hosts_taxon","hosts_gb","hosts_tax_id","hosts_superkingdom","hosts_kingdom","hosts_phylum","hosts_class","hosts_order","hosts_family","hosts_subfamily","hosts_genus","hosts_species"]
#df = df.reindex(columns=column_names)
#df.replace("", "NA", inplace=True)
#df.to_csv(r'contig_host_test.csv', sep=',',index=False)

df1=pd.DataFrame.from_dict(otu_tax,orient='index')
df1.index.name = 'tax_id'
df1.reset_index(inplace=True)
column_names=["tax_id","hosts_taxon","hosts_gb","hosts_tax_id","hosts_superkingdom","hosts_kingdom","hosts_phylum","hosts_class","hosts_order","hosts_family","hosts_subfamily","hosts_genus","hosts_species"]
df1 =df1.reindex(columns=column_names)
df1.replace("", "NA", inplace=True)
df1.to_csv(sys.argv[4], sep=',',index=False)


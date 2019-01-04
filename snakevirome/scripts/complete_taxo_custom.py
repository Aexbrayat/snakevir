import sys, os
from ete3 import NCBITaxa
from collections import Counter
import pandas as pd
import csv
import subprocess

# Taxonomy_file=sys.argv[]
# Taxonomy_full_file=sys.argv[2]
Ref_taxonomy_file=sys.argv[1]
ICTV_taxo_file=sys.argv[2]
Host_taxo_file=sys.argv[3]
output_stat= sys.argv[4]
# output_seq= sys.argv[7]


ref_lineage=[]
ICTV_lineage=[]
host_lineage=[]
host_dict={}
# Taxonomy=open(Taxonomy_file, 'r')
# Taxonomy_full=open(Taxonomy_full_file, 'r')
Ref_taxonomy=open(Ref_taxonomy_file, 'r')
ICTV_taxo=open(ICTV_taxo_file, 'r')
host_taxo=open(Host_taxo_file, 'r')

for lineage in Ref_taxonomy:
	ref_lineage.append(lineage.replace('\n','').split(','))
#
for lineage in host_taxo:
	host_lineage=(lineage.replace('\n','').split(','))
	if host_lineage[0] != "virus_tax_id":
		if host_lineage[0] not in host_dict:
			host_dict[host_lineage[0]]=[host_lineage[1:]]
		else:
			host_dict[host_lineage[0]].append(host_lineage[1:])
for lineage in ICTV_taxo:
	ICTV_lineage.append(lineage.replace('\n','').split(','))

liste_species=[['tax_id','superkingdom','class','clade','order','family','subfamily','genus','species','no rank','taxon','nucleic_acid','hosts_gb','hosts_tax_id','hosts_superkingdom','hosts_class','hosts_order','hosts_family','hosts_subfamily','hosts_genus','hosts_species']]
for lineage in sys.stdin:
	rank=lineage.replace('\n','').replace('\r','').split("	")
	ranks=rank[0:2]+["<not_present>"]+rank[2:]+["<not_present>"]+["<not_present>"]+["<not_present>"]+["<not_present>"]+["<not_present>"]+["<not_present>"]+["<not_present>"]+["<not_present>"]+["<not_present>"]+["<not_present>"]+["<not_present>"]+["<not_present>"]
	liste_species.append(ranks)

seqid=[]

for ranks in liste_species:
	index=liste_species.index(ranks)
	if ranks[1]=="Viruses":
		if ranks[7]!="<not_present>" and ranks[5]=="<not_present>":
			for ref_rank in ref_lineage:
				if ranks[7]==ref_rank[3]:
					liste_species[index][5]=ref_rank[2]
					liste_species[index][4]=ref_rank[1]
					liste_species[index][3]=ref_rank[0]
					liste_species[index][10]=ref_rank[4]
					liste_species[index][11]=ref_rank[5]
		if ranks[5]!="<not_present>" and ranks[4]=="<not_present>":
			for ref_rank in ref_lineage:
				if ranks[5]==ref_rank[2]:
					liste_species[index][4]=ref_rank[1]
					liste_species[index][3]=ref_rank[0]
					liste_species[index][10]=ref_rank[4]
					liste_species[index][11]=ref_rank[5]
		if ranks[5]!="<not_present>" and liste_species[index][10]=="<not_present>":
			liste_species[index][10]=ranks[5]
		if liste_species[index][10]=="<not_present>" and ranks[6]!="<not_present>" :
			liste_species[index][10]=ranks[6]
		if liste_species[index][10]=="<not_present>" and ranks[4]!="<not_present>" :
			liste_species[index][10]=rank[4]
		if liste_species[index][11]=="<not_present>":
			if ranks[8]!="<not_present>":
				for ictv_rank in ICTV_lineage:
					if ranks[7]==ictv_rank[2]:
						liste_species[index][11]=ictv_rank[4]
		if liste_species[index][10]=="<not_present>":
			if "ssRNA negative-strand viruses" in liste_species[index][11]:
				liste_species[index][10]="ssRNA(-)"
			if "ssRNA positive-strand viruses" in liste_species[index][11]:
				liste_species[index][10]="ssRNA(+)"
			if "dsRNA viruses" in liste_species[index][11]:
				liste_species[index][10]="dsRNA"
			if "RNA satellites" in liste_species[index][11]:
				liste_species[index][10]="RNA satellites"
			if "ssDNA viruses" in liste_species[index][11]:
				liste_species[index][10]="ssDNA(+/-)"
			if "dsDNA viruses" in liste_species[index][11]:
				liste_species[index][10]="dsDNA"
			if "unclassified ssRNA viruses" in liste_species[index][11]:
				liste_species[index][10]="ssRNA"
			if "unclassified RNA viruses" in liste_species[index][11]:
				liste_species[index][10]="RNA"
		if liste_species[index][0] not in seqid :
			seqid.append(liste_species[index][0])
#####---------------------------------------------------OK
host_all=[]
ncbi = NCBITaxa()
liste_host=[]
stat_by_seq=[]
n=0
desired_ranks = ['superkingdom', 'class', 'order', 'family', 'subfamily','genus','species']
for ids in seqid:
	n+=1
	ranks2lineage={}
	liste=[]
	fieldnames=[]
	fieldnames1=[]
	host_list=[]
	host_id="<not_present>"
	command = 'esearch -db Nucleotide -query '+ str(ids) +'| efetch -format gb | grep -oP "/host=\K.*"'
	proc = subprocess.Popen(command,stdout=subprocess.PIPE,shell=True)
	(out, err) = proc.communicate()
	host_gb=out.decode("utf-8").rstrip().replace('"','').replace("[u'","").replace("u'","")
	cleaned_host_gb=host_gb.split(' (', 1)[0]
	if "mosquito" in str(cleaned_host_gb.encode('utf8')) or "mosquitoes" in  str(cleaned_host_gb.encode('utf8')):
		cleaned_host_gb="Culicoidea"
	if cleaned_host_gb!='':
		liste.extend([ids , cleaned_host_gb])
		liste_host.append(liste)
		name2taxid = ncbi.get_name_translator([cleaned_host_gb])
		fieldnames.append(host_gb)
		if cleaned_host_gb in name2taxid:
			host_id=int(name2taxid[cleaned_host_gb][0])
			lineage = ncbi.get_lineage(host_id)
			fieldnames.append(host_id)
	if cleaned_host_gb!='' and cleaned_host_gb in name2taxid:
		lineage2name = ncbi.translate_to_names(lineage)
		lineage2ranks=ncbi.get_rank(lineage)
		Dic_lineage2name=dict(zip(lineage, lineage2name))
		for (taxid,rank) in lineage2ranks.items():
			if rank in desired_ranks:
				if not rank in ranks2lineage:
					ranks2lineage[rank]=[taxid]
				else:
					ranks2lineage[rank].append(taxid)
		ranks2names={}
		for rank in ranks2lineage:
			ranks2names[rank]=""
			for i in range(len(ranks2lineage[rank])):
				ranks2names[rank]+=Dic_lineage2name[ranks2lineage[rank][i]]+";"
			ranks2names[rank]=ranks2names[rank][:-1]
		ranks2names['host_id']=host_id
		for rank in desired_ranks:
			fieldnames1.append(ranks2names.get(rank.encode("utf-8"), '<not_present>'))
	host_list=[ids]+fieldnames+fieldnames1
	host_all.append(host_list)

dict_host={}
for seq in liste_species:
	index=liste_species.index(seq)
	if seq[1].isdigit():
		for host in host_all:
			if len(host)>=2 and host[0]==seq[0]:
				for idx_host, val_host in enumerate(host):
					if  isinstance(val_host,str):
						val_host=val_host.encode("utf-8")
					if seq[0] not in dict_host:
						dict_host[seq[0]]={idx_host:[val_host]}
					elif idx_host not in dict_host[seq[0]]:
						dict_host[seq[0]].update({idx_host:[(val_host)]})
					elif  val_host!='<not_present>' and val_host not in dict_host[seq[0]][idx_host]:
						dict_host[seq[0]][idx_host].append(val_host)
			if seq[1] in host_dict:
				host_id=""
				superkingdom=[]
				clas=[]
				order=[]
				family=[]
				subfamily=[]
				genus=[]
				species=[]
				for hosts in host_dict[seq[1]]:
					if host_id=="" and hosts[0]!='':
						host_id+=str(hosts[0]).replace("[","").replace("]","")
					else:
						host_id+="|"+str(hosts[0])
					if hosts[1] not in superkingdom:
						superkingdom.append(hosts[1])
					if hosts[2] not in clas:
						clas.append(hosts[2])
					if hosts[3] not in order:
						order.append(hosts[3])
					if hosts[4] not in family:
						family.append(hosts[4])
					if hosts[5] not in subfamily:
						subfamily.append(hosts[5])
					if hosts[6] not in genus:
						genus.append(hosts[6])
					if hosts[7] not in species:
						species.append(hosts[7])
				liste_species[index][14]=(str(host_id).replace("]","").replace("[","").replace("'","").replace(", ","|"))
				liste_species[index][15]=(str(superkingdom).replace("]","").replace("[","").replace("'","").replace(", ","|"))
				liste_species[index][16]=(str(clas).replace("]","").replace("[","").replace("'","").replace(", ","|"))
				liste_species[index][17]=(str(order).replace("]","").replace("[","").replace("'","").replace(", ","|"))
				liste_species[index][18]=(str(family).replace("]","").replace("[","").replace("'","").replace(", ","|"))
				liste_species[index][19]=(str(subfamily).replace("]","").replace("[","").replace("'","").replace(", ","|"))
				liste_species[index][20]=(str(genus).replace("]","").replace("[","").replace("'","").replace(", ","|"))
				liste_species[index][21]=(str(species).replace("]","").replace("[","").replace("'","").replace(", ","|"))

			else:
				if str(liste_species[index][0]) in dict_host and len(dict_host[str(liste_species[index][0])])>2:
					liste_species[index][14]=str(dict_host[str(liste_species[index][0])][2]).replace("]","").replace("[","").replace("'","").replace(", ","|")
					liste_species[index][15]=str(dict_host[str(liste_species[index][0])][3]).replace("]","").replace("[","").replace("'","").replace(", ","|")
					liste_species[index][16]=str(dict_host[str(liste_species[index][0])][4]).replace("]","").replace("[","").replace("'","").replace(", ","|")
					liste_species[index][17]=str(dict_host[str(liste_species[index][0])][5]).replace("]","").replace("[","").replace("'","").replace(", ","|")
					liste_species[index][18]=str(dict_host[str(liste_species[index][0])][6]).replace("]","").replace("[","").replace("'","").replace(", ","|")
					liste_species[index][19]=str(dict_host[str(liste_species[index][0])][7]).replace("]","").replace("[","").replace("'","").replace(", ","|")
					liste_species[index][20]=str(dict_host[str(liste_species[index][0])][8]).replace("]","").replace("[","").replace("'","").replace(", ","|")
					liste_species[index][21]=str(dict_host[str(liste_species[index][0])][9]).replace("]","").replace("[","").replace("'","").replace(", ","|")
				elif str(liste_species[index][1]) in dict_host and len(dict_host[str(liste_species[index][1])])>=2:
					liste_species[index][13]=str(dict_host[str(liste_species[index][0])][1]).replace("]","").replace("[","").replace("'","").replace(", ","|")
				else:
					continue
		if str(liste_species[index][0]) in dict_host:
			liste_species[index][13]=str(dict_host[str(liste_species[index][0])][1]).replace("]","").replace("[","").replace("'","").replace(", ","|")

list_vir_ids=[]
dic_host_ids={}
to_write=[]
for listes in liste_species:
	if listes[1]=="tax_id":
		to_write.append(listes[0:13]+["host_taxon"]+listes[13:])
	else:
		to_write.append(listes[0:13]+["<not_present>"]+listes[13:])
	if "|" in listes[14] :
		dic_host_ids[listes[0]]=(listes[14].split('|'))
	else :
		dic_host_ids[listes[0]]=([listes[14]])

ranks2lineage={}
ncbi = NCBITaxa()
dic={}
test=[]
for key in dic_host_ids:
	for ids in dic_host_ids[key]:
		if ids != "<not_present>" and ids.isdigit():
			lineage=ncbi.get_lineage(ids)
			lineage2name = ncbi.translate_to_names(lineage)
			lineage2ranks=ncbi.get_rank(lineage)
			Dic_lineage2name=dict(zip(lineage, lineage2name))
			for (taxid,rank) in lineage2ranks.items():
				if taxid in Dic_lineage2name:
					if key not in dic:
						dic[key]={}
					if rank not in dic[key] :
						dic[key][rank]=[Dic_lineage2name[taxid]]
					elif Dic_lineage2name[taxid] not in dic[key][rank]:
						dic[key][rank].append(Dic_lineage2name[taxid])

	if key in dic:
		if 'kingdom' in dic[key] and len(dic[key]['kingdom'])==1 and dic[key]['kingdom'][0]!='<not_present>' and dic[key]['kingdom'][0]!='Metazoa':
			for line in to_write :
				if key in line:
					line[13]=str(dic[key]['kingdom'][0])

		if 'superkingdom' in dic[key] and len(dic[key]['superkingdom'])==1 and dic[key]['superkingdom'][0]!="Eukaryota":
			for line in to_write:
				if key in line:
					line[13]=str(dic[key]['superkingdom'][0])
		if  'order' in dic[key] and 'class' in dic[key] and len(dic[key]['order'])==1 and len(dic[key]['class'])==1:
			if dic[key]['order'][0]=='Diptera':
				if 'family' in dic[key] and dic[key]['family'][0]=='Culicidae':
					for line in to_write:
						if key in line:
							line[13]=str(dic[key]['family'][0])
				else:
					for line in to_write:
						if key in line:
							line[13]=str(dic[key]['order'][0])
		for line in to_write:
			if key in line and  line[13]=='<not_present>' and 'class' in dic[key] and len(dic[key]['class'])==1:
				if dic[key]['class'][0]=='Insecta':
							line[13]=str(dic[key]['phylum'][0])
				elif dic[key]['class'][0]!='Mammalia':
					if 'Vertebrata' in dic[key]['no rank']:
								line[13]='Vertebrata'
				elif dic[key]['class'][0]=='Mammalia':
							line[13]=str(dic[key]['class'][0])
			if key in line and  'phylum' in dic[key] and  line[11]=='<not_present>' and dic[key]['phylum'][0]=='Arthropoda':
						line[13]='Arthropoda'
			if key in line and  line[13]=='<not_present>':
				if len(line[13].split('|'))==1:
					line[13]='Other'
with open(output_stat, "w") as f:
    writer = csv.writer(f)
    writer.writerows(to_write)
# with open(output_seq, "w") as g:
#     writer = csv.writer(g)
#     writer.writerows(stat_by_seq)
# # f.close
# # g.close

import sys, os
from collections import Counter
import pandas as pd
import csv
import subprocess
import pickle

with open(sys.argv[1], 'rb') as handle:
	taxo_dict = pickle.load(handle)
with open(sys.argv[2], 'rb') as handle:
	taxo_custom = pickle.load(handle)


Taxonomy_file=sys.argv[3]
output_stat=sys.argv[4]

Taxonomy=open(Taxonomy_file, 'r')



liste_species=[]
for lineage in Taxonomy:
	rank=lineage.replace('\n','').replace('\r','').split(",")
	
	if rank[0]=="tax_id":
		ranks=rank[0:1]+["Nucleic_acid"]+["genome"]+["cluster"]+rank[1:2]+rank[3:]+["how"]
		print(ranks)
	elif rank[2]=="Viruses":
		ranks=rank[0:1]+["NA"]+["NA"]+["NA"]+rank[1:2]+rank[3:]+["NA"]
	else:
		continue
	liste_species.append(ranks)

for ranks in liste_species:
	index=liste_species.index(ranks)

def complet_ictv(liste_species):
	for ranks in liste_species:
		index=liste_species.index(ranks)

		if liste_species[index][11] in taxo_dict['species'] and liste_species[index][10]=='NA':
			liste_species[index][10]=taxo_dict['species'][liste_species[index][11]]

		if liste_species[index][10] in taxo_dict['genus'] and liste_species[index][9]=='NA':
			liste_species[index][9]=taxo_dict['genus'][liste_species[index][10]]
			
		if liste_species[index][9]in taxo_dict['family'] and liste_species[index][8]=='NA':
			liste_species[index][8]=taxo_dict['family'][liste_species[index][9]]

			
		if liste_species[index][8] in taxo_dict['order'] and liste_species[index][7]=='NA':
			liste_species[index][7]=taxo_dict['order'][liste_species[index][8]]

			
		if liste_species[index][7] in taxo_dict['clas'] and liste_species[index][6]=='NA':
			liste_species[index][6]=taxo_dict['clas'][liste_species[index][7]]

			
		if liste_species[index][6] in taxo_dict['phylum'] and liste_species[index][5]=='NA':
			liste_species[index][5]=taxo_dict['phylum'][liste_species[index][6]]

			
		if liste_species[index][5] in taxo_dict['kingdom'] and liste_species[index][4]=='NA':
			liste_species[index][4]=taxo_dict['kingdom'][liste_species[index][5]]

		if liste_species[index][9] in taxo_dict['fam2gen'] and liste_species[index][4]=='NA':
			liste_species[index][4]=taxo_dict['fam2gen'][liste_species[index][9]]
	
complet_ictv(liste_species)
print("ok ictv")
n=0
for ranks in liste_species:
	n+=1
	print(n,"/",len(liste_species))
	index=liste_species.index(ranks)
	if liste_species[index][10] not in taxo_dict['species'] and liste_species[index][11] in taxo_dict['taxshi']:
		if liste_species[index][10]=="NA" and  taxo_dict['taxshi'][liste_species[index][11]][6]!="NA":
			liste_species[index][10]=taxo_dict['taxshi'][liste_species[index][11]][6]
			liste_species[index][13]+=" | shi_genus "
		if liste_species[index][9]=="NA" and  taxo_dict['taxshi'][liste_species[index][11]][5]!="NA":
			liste_species[index][9]=taxo_dict['taxshi'][liste_species[index][11]][5]
			liste_species[index][13]+=" | shi_fam "
		if liste_species[index][8]=="NA" and  taxo_dict['taxshi'][liste_species[index][11]][4]!="NA":
			liste_species[index][8]=taxo_dict['taxshi'][liste_species[index][11]][4]
			liste_species[index][13]+=" | shi_order "
		if liste_species[index][7]=="NA" and  taxo_dict['taxshi'][liste_species[index][11]][3]!="NA":
			liste_species[index][7]=taxo_dict['taxshi'][liste_species[index][11]][3]
			liste_species[index][13]+=" | shi_class "
		if liste_species[index][6]=="NA" and  taxo_dict['taxshi'][liste_species[index][11]][2]!="NA":
			liste_species[index][6]=taxo_dict['taxshi'][liste_species[index][11]][2]
			liste_species[index][13]+=" | shi-phyl "
		if liste_species[index][5]=="NA" and  taxo_dict['taxshi'][liste_species[index][11]][1]!="NA":
			liste_species[index][5]=taxo_dict['taxshi'][liste_species[index][11]][1]
			liste_species[index][13]+=" | shi_king "
		if liste_species[index][4]=="NA" and  taxo_dict['taxshi'][liste_species[index][11]][0]!="NA":
			liste_species[index][4]=taxo_dict['taxshi'][liste_species[index][11]][0]
			liste_species[index][13]+=" | shi_clade "
	complet_ictv(liste_species)
	
	if liste_species[index][11] not in taxo_dict['species'] and liste_species[index][11] in taxo_custom['taxcustom']:
		
		if liste_species[index][10]=="NA" and  taxo_custom['taxcustom'][liste_species[index][11]][10]!="NA":
			liste_species[index][10]=taxo_custom['taxcustom'][liste_species[index][11]][10]
			liste_species[index][13]+=" | custom_genus "
			
		if liste_species[index][9]=="NA" and  taxo_custom['taxcustom'][liste_species[index][11]][8]!="NA":
			liste_species[index][9]=taxo_custom['taxcustom'][liste_species[index][11]][8]
			liste_species[index][13]+=" | custom_fam "
			
		if liste_species[index][8]=="NA" and  taxo_custom['taxcustom'][liste_species[index][11]][7]!="NA":
			liste_species[index][8]=taxo_custom['taxcustom'][liste_species[index][11]][7]
			liste_species[index][13]+=" | custom_order "
			
		if liste_species[index][7]=="NA" and  taxo_custom['taxcustom'][liste_species[index][11]][6]!="NA":
			liste_species[index][7]=taxo_custom['taxcustom'][liste_species[index][11]][6]
			liste_species[index][13]+=" | custom_class "
			
		if liste_species[index][6]=="NA" and  taxo_custom['taxcustom'][liste_species[index][11]][5]!="NA":
			liste_species[index][6]=taxo_custom['taxcustom'][liste_species[index][11]][5]
			liste_species[index][13]+=" | custom_phyl "
			
		if liste_species[index][5]=="NA" and  taxo_custom['taxcustom'][liste_species[index][11]][4]!="NA":
			liste_species[index][5]=taxo_custom['taxcustom'][liste_species[index][11]][4]
			liste_species[index][13]+=" | custom_king "
			
		if liste_species[index][4]=="NA" and  taxo_custom['taxcustom'][liste_species[index][11]][3]!="NA":
			liste_species[index][4]=taxo_custom['taxcustom'][liste_species[index][11]][3]
			liste_species[index][13]+=" | custom_clade "
		if liste_species[index][3]=="NA":
			liste_species[index][3]=taxo_custom['taxcustom'][liste_species[index][11]][2]

			liste_species[index][13]+=" | custom_clust "
		
	complet_ictv(liste_species)
		
	if liste_species[index][11] in taxo_dict['clustictv'] and liste_species[index][3]=='NA':
		liste_species[index][3]=taxo_dict['clustictv'][liste_species[index][11]]

	
	if liste_species[index][3]=='NA' and liste_species[index][9]!="NA":
		liste_species[index][3]=(liste_species[index][9].replace('viridae', ''))

	if liste_species[index][11] in taxo_dict['clustshi'] and liste_species[index][3]=='NA':
		liste_species[index][3]=taxo_dict['clustshi'][liste_species[index][11]]
		liste_species[index][13]+=" | clustshi "
		
	if liste_species[index][3]=='NA' and liste_species[index][8]!="NA":
		liste_species[index][3]=(liste_species[index][8].replace('virales', '')+"_unclass")
		
	if liste_species[index][1]=='NA' and liste_species[index][7] in taxo_dict['clas2gent']:
		liste_species[index][1]=taxo_dict['clas2gent'][liste_species[index][7]]
	if liste_species[index][1]=='NA' and liste_species[index][8] in taxo_dict['ord2gent']:
		liste_species[index][1]=taxo_dict['ord2gent'][liste_species[index][8]]
	if liste_species[index][1]=='NA' and liste_species[index][9] in taxo_dict['fam2gent']:
		liste_species[index][1]=taxo_dict['fam2gent'][liste_species[index][9]]

	if liste_species[index][1]=="NA":
		if "ssRNA negative-strand viruses" in liste_species[index][12] or "ssRNA negative-strand viruses" in   liste_species[index][4]:
			liste_species[index][1]="ssRNA(-)"
		if "ssRNA positive-strand viruses" in liste_species[index][12]or "ssRNA positive-strand viruses"  in   liste_species[index][4]:
			liste_species[index][1]="ssRNA(+)"
		if "dsRNA" in liste_species[index][12]or "dsRNA" in   liste_species[index][4]:
			liste_species[index][1]="dsRNA"
		if "RNA satellites" in liste_species[index][12] or "RNA satellites" in liste_species[index][4]:
			liste_species[index][1]="RNA satellites"
		if "ssDNA" in liste_species[index][12] or  "ssDNA"  in liste_species[index][4]:
			liste_species[index][1]="ssDNA(+/-)"
		if "dsDNA" in liste_species[index][12]  or  "dsDNA" in liste_species[index][4]:
			liste_species[index][1]="dsDNA"
		if "unclassified ssRNA" in liste_species[index][12] or "unclassified ssRNA"in liste_species[index][4]:
			liste_species[index][1]="ssRNA"
		if "unclassified RNA" in liste_species[index][12] or "unclassified RNA" in liste_species[index][4]:
	 		 liste_species[index][1]="RNA"
	if liste_species[index][2]=="NA":
		if 'DNA' in liste_species[index][1] :
			liste_species[index][2]="DNA"
		if 'RNA' in liste_species[index][1] :
			liste_species[index][2]="RNA"
	if liste_species[index][2]=="NA":
		if 'DNA' in liste_species[index][11] :
			liste_species[index][2]="DNA"
		if 'RNA' in liste_species[index][11] :
			liste_species[index][2]="RNA"
	
	if liste_species[index][2]=="NA":
		if "Riboviria" in liste_species[index][4]:
			liste_species[index][2]="RNA"
		if "Varidnaviria" in liste_species[index][4]:
			liste_species[index][2]="DNA"
		if "Duplodnaviria" in liste_species[index][4] or "Monodnaviria"in liste_species[index][4] :
			liste_species[index][2]="DNA"
			
	if liste_species[index][4]=="NA":
		if "RNA" in liste_species[index][2]:
			liste_species[index][4]="Riboviria"
		if "DNA" in liste_species[index][2]:
			liste_species[index][4]="Varidnaviria"
		if "DNA" in liste_species[index][4]:
			liste_species[index][4]="Duplodnaviria"
			
	if liste_species[index][6]=='NA' and liste_species[index][1]!='NA':
		liste_species[index][6]=liste_species[index][2]

print(liste_species)
with open(output_stat, "w") as f:
    writer = csv.writer(f)
    writer.writerows(liste_species)


import sys, os
from collections import Counter
import pandas as pd
import csv

#Take 3 input Blast resum , Mat files , lineage
#Output two files Coverage  of contigs for each samples and coverage for lineages.

#Files with the blast infos for each seq/contigs. Output of the rule Join_seq_acc_taxo
Seq_hit_info_file= sys.argv[1]
with open(Seq_hit_info_file,'r') as f:
	seq_hit_info = f.readlines()[1:]
#seq_hit_info = open(Seq_hit_info_file, 'r')

# Path to the directory CountsMapping_all/. Outputs of the rule Quantify_contigs_coverage
Path_file_mat= sys.argv[2]
file_mat= os.listdir(Path_file_mat)

Path_file_mat_ext=sys.argv[3]

# File with the taxonomic lineage for each tax_ids. Output of the rule get_lineage_from_taxids
Taxonomy_file=sys.argv[4]
Taxonomy=open(Taxonomy_file, 'r')


output_full= sys.argv[5]


Dic_count={}
for files in file_mat:
	if files.endswith(Path_file_mat_ext):
		f=open(Path_file_mat+files, 'r')
		files_name=files.split(Path_file_mat_ext)[0]
		for line in f :
			line=line.replace('\n','').split()
			print(line)
			if files_name not in Dic_count and files_name!="" :
				Dic_count[files_name]={}
				Dic_count[files_name].update({line[1]:line[0]})
			else:
				Dic_count[files_name].update({line[1]:line[0]})
data_frame = pd.DataFrame(Dic_count).fillna(0)
del Dic_count
del line
f.close()



Dic_count_by_samples={}
data_frame_to_str=data_frame.to_csv()
list_count=data_frame_to_str.split("\n")
samples=list_count[0].split(",")

for line in list_count[1:]:
	liste_count=line.split(',')
	print(liste_count)
	Dic_count_by_samples[liste_count[0]]=[]
	for count in liste_count[1:]:
		Dic_count_by_samples[liste_count[0]].append(int(count))
del list_count
del data_frame_to_str


### Produce a dic with Tax_ids as key and list of contig/seq corresponding as value
Dic_best_hit={}
for hits in seq_hit_info:
	list_hit=hits.replace('\n','').split(',')
	print(list_hit)
	list_hit[0]=list_hit[0].split('/', 1)[0]
	if list_hit[-1] not in Dic_best_hit :
		Dic_best_hit[list_hit[-1]]={list_hit[0]:list_hit[1:-1]}
	else :
		Dic_best_hit[list_hit[-1]].update({list_hit[0]:list_hit[1:-1]})

print(Dic_best_hit)
Dic_count_taxids={}
Sum_count_per_taxids={}
for tax_ids in Dic_best_hit:
	l=[]
	Dic_count_taxids[tax_ids]=[]
	for ids in Dic_best_hit[tax_ids]:
		ids=ids.split('/', 1)[0] #Remove the /1 or /2 from the ids of sequences
		try:
			Dic_count_taxids[tax_ids].append(Dic_count_by_samples[ids])
		except KeyError:
			continue
	for listes in Dic_count_taxids[tax_ids]:
		l.append(listes)
		Sum_count_per_taxids[tax_ids]=list(map(sum, zip(*l)))

#tax_hits_infos={}
#for tax_ids in Dic_best_hit:
#	tax_hits_infos[tax_ids]=[]
#	num_contigs=0
#	total_query_length_cont=0
#	total_subject_length=0
#	total_align_length_cont=0
#	ratio_length=0
#	percentage_id=0
#	percent_query=0
#	percent_subject=0
#	match_length=0
#	list_subject=[]
#	number_sigleton=0
#	total_align_length_sin=0
#	percentage_id_sigleton=0
#	contig_id=[]
#	average_p_id=0
#	for seqs in Dic_best_hit[tax_ids]:
#		if seqs.startswith("k") or seqs.startswith("C") :
#			Dic_best_hit[tax_ids][seqs]
#			num_contigs+=1
#			total_align_length_cont+=int(Dic_best_hit[tax_ids][seqs][3])
#			total_query_length_cont+=int(Dic_best_hit[tax_ids][seqs][1])
#			total_subject_length+=int(Dic_best_hit[tax_ids][seqs][2])
#			percent_query+=float(Dic_best_hit[tax_ids][seqs][8])*((int(Dic_best_hit[tax_ids][seqs][1])))
#			percent_subject+=((float(Dic_best_hit[tax_ids][seqs][3])*100)/int(Dic_best_hit[tax_ids][seqs][2]))*((int(Dic_best_hit[tax_ids][seqs][2])))
#			percentage_id+=float(Dic_best_hit[tax_ids][seqs][9])*(float(Dic_best_hit[tax_ids][seqs][3]))
#			match_length+=(float(Dic_best_hit[tax_ids][seqs][3]))
#			average_p_id+=(float(Dic_best_hit[tax_ids][seqs][9]))
#			contig_id.append(float(Dic_best_hit[tax_ids][seqs][9]))
#			if Dic_best_hit[tax_ids][seqs][0] not in list_subject:
#				list_subject.append(Dic_best_hit[tax_ids][seqs][0])
#			else:
#				continue
#		else:
#			number_sigleton+=1
#			total_align_length_sin+=int(Dic_best_hit[tax_ids][seqs][3])
#			if number_sigleton!=0:
#				percentage_id_sigleton+=float(Dic_best_hit[tax_ids][seqs][9])*((int(Dic_best_hit[tax_ids][seqs][3])))
#			else:
#				percentage_id_sigleton=0
#	if num_contigs !=0 and number_sigleton!=0 :
#		tax_hits_infos[tax_ids].extend(( "%.2f" % float(match_length/num_contigs), num_contigs, number_sigleton, "%.2f" % float(total_query_length_cont/num_contigs), total_query_length_cont, "%.2f" % ((percent_query/(int(total_query_length_cont)))), len(list_subject), "%.2f" % ((percent_subject/(int(total_subject_length)))), "%.2f" % ((percentage_id/(int(total_align_length_cont)))), min(contig_id) , max(contig_id) , "%.2f" % (average_p_id/num_contigs),
#		"%.2f" % ((percentage_id_sigleton/(int(total_align_length_sin)))/number_sigleton)))
#	if num_contigs !=0 and number_sigleton==0 :
#		tax_hits_infos[tax_ids].extend(("%.2f" % float(match_length/num_contigs), num_contigs, number_sigleton, "%.2f" % float(total_query_length_cont/num_contigs), total_query_length_cont, "%.2f" % ((percent_query/(int(total_query_length_cont)))/num_contigs), len(list_subject), "%.2f" % ((percent_subject/(int(total_subject_length)))), "%.2f" % ((percentage_id/(int(total_align_length_cont)))),min(contig_id) , max(contig_id) , "%.2f" % (average_p_id/num_contigs), "-"))
#	elif num_contigs ==0 and number_sigleton !=0 :
#		tax_hits_infos[tax_ids].extend(( "-", "-" , number_sigleton, "-", "-","-", "-", "-", "-","-","-","-", "%.2f" % ((percentage_id_sigleton/(int(total_align_length_sin))))))
#	tax_hits_infos[tax_ids]=tax_hits_infos[tax_ids]+Sum_count_per_taxids[tax_ids][:]

ranks=Taxonomy.readline().replace('\n','').split("|")
#stat_info=["Avg_match_length","n_Contigs","n_Singleton","Avg_contigs_length","Sum_contigs_length","Avg_Coverage_contigs","n_subject","Avg_Coverage_subject","Pond_Avg_contigs_p_id","min_contigs_p_id","max_contig_p_id", "avg_contig_p_id","Avg_singleton_p_id" ]
#lineage_coverage=[ranks+stat_info+samples[1:]]

stat_info_brut=["qseqid"]
lineage_brut=[stat_info_brut+samples[1:]]
dic_lineage={}
#for lineage in Taxonomy:
#	list_lineage=lineage.replace('\n','').split('|')
#	lineage_coverage.append(list_lineage + tax_hits_infos[list_lineage[0]])
#	dic_lineage[list_lineage[0]]=list_lineage

list_table=[]
DIc_test={}
for tax_id in Dic_best_hit:
	for key in Dic_best_hit[tax_id].keys():
		key=key.split('/', 1)[0]
		try:
			lineage_brut.append([key]+list_table + Dic_count_by_samples[key])
		except KeyError:
			continue
#with open(output_stat, "w") as f:
#    writer = csv.writer(f)
#    writer.writerows(lineage_coverage)
with open(output_full, "w") as f1:
    writer = csv.writer(f1)
    
    writer.writerows(lineage_brut)
Taxonomy.close()

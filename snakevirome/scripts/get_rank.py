# #!/bin/python

import csv
from ete3 import NCBITaxa
import sys
from collections import Counter

ncbi = NCBITaxa()

def get_desired_ranks(taxids, desired_ranks):
	lineage = ncbi.get_lineage(taxids)
	lineage2ranks = ncbi.get_rank(lineage)
	ranks2lineage = {}
	for (taxid,rank) in lineage2ranks.items():
		if not rank in ranks2lineage:
			ranks2lineage[rank]=[taxid]
		else:
			ranks2lineage[rank].append(taxid)
	lineage2name = ncbi.translate_to_names(lineage)
	Dic_lineage2name=dict(zip(lineage, lineage2name))

	ranks2names={}
	for rank in ranks2lineage:
		ranks2names[rank]=""
		for i in range(len(ranks2lineage[rank])):
			ranks2names[rank]+=Dic_lineage2name[ranks2lineage[rank][i]]+";"
		ranks2names[rank]=ranks2names[rank][:-1]
	ranks2names['tax_id']=taxids

	return {'{}'.format(rank): ranks2names.get(rank, '<not_present>') for rank in desired_ranks}

def main(taxids, desired_ranks, path):
    with open(path, 'w') as csvfile:
        fieldnames = ['{}'.format(rank) for rank in desired_ranks]
        writer = csv.DictWriter(csvfile, delimiter='|', fieldnames=fieldnames)
        writer.writeheader()
        for taxid in taxids:
            taxid=taxid.replace('\n','')
            writer.writerow(get_desired_ranks(taxid, desired_ranks))
    csvfile.close()


if __name__ == '__main__':
    taxids=sys.stdin.read().split(',')
    desired_ranks = ['tax_id', 'superkingdom', 'class', 'order', 'family', 'subfamily','genus','species', 'no rank']
    path = sys.argv[1]
    main(taxids, desired_ranks, path)

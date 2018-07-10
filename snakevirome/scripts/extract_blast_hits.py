# #!/bin/python

import sys, os
from Bio import SeqIO
blast_file = sys.argv[1]
fasta_file = sys.argv[2]
output_file = sys.argv[3]



id_blast=[]
with open(blast_file) as f:
    for line in f:
        id_blast.append(line.split()[0])
    id_blast=set(id_blast)
    id_blast=list(id_blast)

fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
with open(output_file, "w") as f:
    for seq in fasta_sequences:
        if seq.id in id_blast:
            SeqIO.write([seq], f, "fasta")

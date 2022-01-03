import sys, os
from collections import Counter
import pandas as pd
import csv
from os.path import exists
import glob
import numpy as np

Path_file_mat= sys.argv[1]
file_mat= os.listdir(Path_file_mat)
output_full= sys.argv[2]

count_table = pd.DataFrame({'qseqid' : []})
n=0
for files in file_mat:
	print(n)
	n+=1
	if "contig" in files:
		df=pd.read_csv(Path_file_mat+files, sep=' ').drop_duplicates()
		count_table=count_table.merge(df, on='qseqid', how='outer').fillna(0)
count_table.iloc[: , 1:] = count_table.iloc[: , 1:].astype(int)

count_table.to_csv(sys.argv[2],index=False)

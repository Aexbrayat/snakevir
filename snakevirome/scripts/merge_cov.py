import sys
import csv
import glob
import pandas as pd
import gc
import os
columns = []
data = {}
ids = set()
df=pd.DataFrame()
for subdir, dirs, files in os.walk(sys.argv[1]):
	for filename in files:
		filepath = subdir + os.sep + filename
		if filepath.endswith(".cov"):
			key = filename.rstrip("_on_nova_bat.bam.cov")
			with open(filepath) as tsvfile:
				f = pd.read_csv(tsvfile, delimiter="\t", names=["contig", "nuc_pos", key])
				if df.empty:
					df=f
				else:
					df=df.merge(f, how='outer', on=['contig','nuc_pos'])
				del f
				gc.collect()
df = df.fillna(0).drop_duplicates()

df.to_csv(sys.argv[2], index=False)

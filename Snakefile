import re
import sys
from os.path import join

def message(mes):
	sys.stderr.write("|---- " + mes + "\n")

def errormes(mes):
	sys.stderr.write("| ERROR ----" + mes + "\n")

configfile : "config.yaml"
workdir: "/scif/data"

datadir = config["fastq"]
scriptdir = config["Scripts"]
filter_2=config["filter_2"]
filter_1=config["filter_1"]
use_filter2=config["use_filter2"]
base_nr = config["base_nr"]
base_taxo = config["base_taxo"]

ext_R1=config["ext_R1"]
ext_R2=config["ext_R2"]
ext=config["ext"]

READS, =  glob_wildcards(datadir+"{readfile}"+ext)
SAMPLES, = glob_wildcards(datadir+"{sample}"+ext_R1+ext)
message(str(use_filter2))
NBSAMPLES = len(SAMPLES)
NBREADS = len(READS)
ref_taxo=config["ref_taxo"]
ref_taxo1=config["ref_taxo1"]
ICTV_taxo=config["ICTV_taxo"]
host_taxo=config["host_taxo"]
RUN = config["run"]
message(str(NBSAMPLES)+" samples  will be analysed")
message(str(len(READS))+" fastq files  will be processed")
message("Run name: "+RUN)
if NBREADS != 2*NBSAMPLES:
	errormes("Please provide two reads file per sample")
	sys.exit()

rule all:
	input:
		# "/Coverage/lineage_revise_stats_"+RUN+".csv",
		expand("logs/insert_size/{smp}_insert_size_metrics_"+RUN+".txt", smp=SAMPLES),
		expand("logs/logs_coverage/{smp}_coverage_"+RUN+".txt", smp=SAMPLES),
		expand("logs/logs_contaminent/Stats_contaminent_{smp}.txt", smp=SAMPLES),
		"logs/Summary_log_"+RUN+".csv",
		"logs/logsAssembly/"+RUN+"_Assembly.log",
		"logs/Summary_log_"+RUN+".csv",
		expand("FiltersMapping/unmapped_{smp}_filter2_widows.sam", smp=SAMPLES),

# Remove sequencing adapters based on 5' and 3' sequences (A3 and A5 in the config.yaml file)
rule Remove_sequencing_adapters:
	input:
		R1=datadir+"{smp}"+ ext_R1+ext,
		R2=datadir+"{smp}"+ ext_R2+ext,
	output:
		R1="cutadaptfiles/{smp}_1.trimmed.fastq",
		R2="cutadaptfiles/{smp}_2.trimmed.fastq",
	log:
		"logs/logscutadapt/{smp}_cut1.log",
	params:
		A3 = config["A3"],
		A5 = config["A5"],
	shell:
		"""
		scif --quiet run cutadapt -n 10 -g {params.A5} -a {params.A3} --overlap 15  -o {output.R1} -p {output.R2}  {input.R1} {input.R2} &> {log}
		"""
# Quality trimming based on the --quality-cutoff (low-quality ends) -q parameter to change the cutoff
# and --minimum-length -m option to Discard processed reads that are shorter than the length specified.
rule Quality_trimming:
	input:
		R1="cutadaptfiles/{smp}_1.trimmed.fastq",
		R2="cutadaptfiles/{smp}_2.trimmed.fastq",
	output:
		R1="cutadaptfiles/{smp}_1.clean.fastq",
		R2="cutadaptfiles/{smp}_2.clean.fastq",
	log:
		"logs/logscutadapt/{smp}_cut2.log"
	shell:
		"""
		scif run cutadapt -q 30,30 -m 40 -o {output.R1} -p {output.R2} {input.R1} {input.R2} &> {log}
		"""

# Resynchronize 2 fastq or fastq.gz files (R1 and R2) after they have been trimmed and cleaned.
rule Repair_Pairs:
	input:
		R1="cutadaptfiles/{smp}_1.clean.fastq",
		R2="cutadaptfiles/{smp}_2.clean.fastq",
	output:
		R1="cutadaptfiles/{smp}_1.clean.fastq_pairs_R1.fastq",
		R2="cutadaptfiles/{smp}_2.clean.fastq_pairs_R2.fastq",
		WI="cutadaptfiles/{smp}_1.clean.fastq_singles.fastq"
	params:
		scriptdir+"fastqCombinePairedEnd.py"
	log:
		"logs/logsRepairspairs/{smp}_repair.log"
	shell:
		"""
		python2 {params} {input.R1} {input.R2}  2> log
		"""

# Filter_1-homologous sequences cleaning by mapping (bwa) on ribosomal sequences.
rule Map_On_filter_1:
	input:
		R1="cutadaptfiles/{smp}_fastq_pairs_R1.fastq",
		R2="cutadaptfiles/{smp}_fastq_pairs_R2.fastq",
		WI="cutadaptfiles/{smp}_fastq_widows.fastq"
	output:
		bam_pairs="/scif/data/FiltersMapping/{smp}_filter1_pairs.bam",
		bam_WI="/scif/data/FiltersMapping/{smp}_filter1_widows.bam",
		sam_pairs="FiltersMapping/unmapped_{smp}_filter1_pairs.sam",
		sam_WI="FiltersMapping/unmapped_{smp}_filter1_widows.sam"
	params:
		filter_1 = config["filter_1"],
		filter_1_index = config["filter_1"]+".bwt",
	log:
		pairs="/scif/data/logs/logsMapFilters/{smp}_bwa_pairs_on_filter1.log",
		WI="/scif/data/logs/logsMapFilters/{smp}_bwa_widows_on_filter1.log"
	threads: 5
	shell:
		"""
		if [ -s {input.R1} ] && [ -s {input.R2} ]
		then
			scif --quiet run bwa mem -t {threads} {params.filter_1} {input.R1} {input.R2} 2> {log.pairs} | samtools view -b - | tee {output.bam_pairs} | scif run samtools view -f 0x4 -o {output.sam_pairs}
		else
			echo "{input.R1} or {input.R2} is empty."
			touch {output.bam_pairs}
			touch {output.sam_pairs}
		fi
		if [ -s {input.WI} ]
		then
		scif --quiet run bwa mem -t {threads} {params.filter_1} {input.WI} 2> {log.WI} | samtools view -b - | tee {output.bam_WI} | scif run samtools view -f 0x4 -o {output.sam_WI}
		else
			echo "{input.WI} is empty."
			touch {output.bam_WI}
			touch {output.sam_WI}
		fi
		"""
# Extact samples squences with with no homologous sequences with the filter_1.
rule Extract_Unmapped_filter1_Reads:
	input:
		sam_pairs="FiltersMapping/unmapped_{smp}_filter1_pairs.sam",
		sam_WI="FiltersMapping/unmapped_{smp}_filter1_widows.sam"
	output:
		pairs_R1="FilteredFastq/filtered_filter1_{smp}_R1.fastq",
		pairs_R2="FilteredFastq/filtered_filter1_{smp}_R2.fastq",
		WI="FilteredFastq/filtered_filter1_{smp}_widows.fastq",
	shell:
		"""
		if [ -s {input.sam_pairs} ]
		then
			scif --quiet run picard SamToFastq VALIDATION_STRINGENCY=SILENT I={input.sam_pairs} F={output.pairs_R1} F2={output.pairs_R2}
		else
			echo "{input.sam_pairs} or {input.sam_WI} is empty."
			touch {output.pairs_R1}
			touch {output.pairs_R2}
		fi
		if [ -s {input.sam_WI} ]
		then
			scif --quiet run picard SamToFastq VALIDATION_STRINGENCY=SILENT I={input.sam_WI} F={output.WI}
		else
			echo "{input.sam_WI} is empty."
			touch {output.WI}
		fi


		"""
# filter2-homologous sequences cleaning by mapping (bwa) on ribosomal sequences.
rule Map_On_filter2:
	input:
		R1="FilteredFastq/filtered_filter1_{smp}_R1.fastq",
		R2="FilteredFastq/filtered_filter1_{smp}_R2.fastq",
		WI="FilteredFastq/filtered_filter1_{smp}_widows.fastq"
	output:
		bam_pairs="/scif/data/FiltersMapping/{smp}_filter2_pairs.bam",
		bam_WI="/scif/data/FiltersMapping/{smp}_filter2_widows.bam",
		sam_pairs="FiltersMapping/unmapped_{smp}_filter2_pairs.sam",
		sam_WI="FiltersMapping/unmapped_{smp}_filter2_widows.sam"

	params:
		filter_2 = config["filter_2"],
		filter_2_index = config["filter_2"]+".bwt",
		use_filter2=config["use_filter2"]
	log:
		pairs="/scif/data/logs/logsMapFilters/{smp}_bwa_pairs_on_filter2.log",
		WI="/scif/data/logs/logsMapFilters/{smp}_bwa_widows_on_filter2.log"
	threads: 5
	shell:
		"""
		if [ "{params.use_filter2}"="yes" ]
		then
			if [ -s {input.R1} ] && [ -s {input.R2} ]
			then
				scif --quiet run bwa mem -t {threads} {params.filter_2} {input.R1} {input.R2} 2> {log.pairs} | samtools view -b - | tee {output.bam_pairs} | scif run samtools view -f 0x4 -o {output.sam_pairs}
			else
				echo "{input.R1} or {input.R2} is empty."
				touch {output.bam_pairs}
				touch {output.sam_pairs}
			fi
			if [ -s {input.WI} ]
			then
			scif --quiet run bwa mem -t {threads} {params.filter_2} {input.WI} 2> {log.WI} | samtools view -b - | tee {output.bam_WI} | scif run samtools view -f 0x4 -o {output.sam_WI}
			else
				echo "{input.WI} is empty."
				touch {output.bam_WI}
				touch {output.sam_WI}
			fi
		fi
		"""

# Extact samples squences with with no homologous sequences with filter2.
rule Extract_Unmapped_MapFilters_Reads:
	input:
		sam_pairs="FiltersMapping/unmapped_{smp}_filter2_pairs.sam",
		sam_WI="FiltersMapping/unmapped_{smp}_filter2_widows.sam",
		R1="FilteredFastq/filtered_filter1_{smp}_R1.fastq",
		R2="FilteredFastq/filtered_filter1_{smp}_R2.fastq",
		WI="FilteredFastq/filtered_filter1_{smp}_widows.fastq"
	output:
		pairs_R1="FilteredFastq/filtered_filter2_{smp}_R1.fastq",
		pairs_R2="FilteredFastq/filtered_filter2_{smp}_R2.fastq",
		WI="FilteredFastq/filtered_filter2_{smp}_widows.fastq"

	params:
		use_filter2=config["use_filter2"]

	shell:
		"""
		if [ "{params.use_filter2}"="yes" ]
		then
			if [ -s {input.sam_pairs} ]
			then
				scif --quiet run picard SamToFastq VALIDATION_STRINGENCY=SILENT I={input.sam_pairs} F={output.pairs_R1} F2={output.pairs_R2}
			else
				echo "{input.sam_pairs} or {input.sam_WI} is empty."
				touch {output.pairs_R1}
				touch {output.pairs_R2}
			fi
			if [ -s {input.sam_WI} ]
			then
				scif --quiet run picard SamToFastq VALIDATION_STRINGENCY=SILENT I={input.sam_WI} F={output.WI}
			else
				echo "{input.sam_WI} is empty."
				touch {output.WI}
			fi
		else
			cat {input.R1} > {output.pairs_R1}
			cat {input.R2} > {output.pairs_R2}
			cat {input.WI} > {output.WI}
		fi
		"""

# Extract some stats (logs) from the cleaning steps
rule log_map_conta:
	input:
		filter1_pairs_R1="FilteredFastq/filtered_filter1_{smp}_R1.fastq",
		filter1_pairs_R2="FilteredFastq/filtered_filter1_{smp}_R2.fastq",
		filter1_WI="FilteredFastq/filtered_filter1_{smp}_widows.fastq",
		filter2_pairs_R1="FilteredFastq/filtered_filter2_{smp}_R1.fastq",
		filter2_pairs_R2="FilteredFastq/filtered_filter2_{smp}_R2.fastq",
		filter2_WI="FilteredFastq/filtered_filter2_{smp}_widows.fastq"
	output:
		"logs/logs_contaminent/Stats_contaminent_{smp}.txt",
	shell:
		"""
		awk '{{s++}}END{{print "filter1_pair_R1 : " s/4}}' {input.filter1_pairs_R1} >> {output};
		awk '{{s++}}END{{print "filter1_pairs_R2 : " s/4}}' {input.filter1_pairs_R2} >> {output};
		awk '{{s++}}END{{print "filter1_pair_wi : " s/4}}' {input.filter1_WI} >> {output};
		awk '{{s++}}END{{print "filter2_pair_R1 : " s/4}}' {input.filter2_pairs_R1} >> {output};
		awk '{{s++}}END{{print "filter2_pairs_R2 : " s/4}}' {input.filter2_pairs_R2} >> {output};
		awk '{{s++}}END{{print "filter2_pair_wi : " s/4}}' {input.filter2_WI} >> {output};
		"""

# Merging overlapping pairs
rule Merge_Pairs_With_Flash:
	input:
		pairs_R1="FilteredFastq/filtered_filter2_{smp}_R1.fastq",
		pairs_R2="FilteredFastq/filtered_filter2_{smp}_R2.fastq",
	params:
		prefix="FilteredFastq/{smp}",
	output:
		ext="FilteredFastq/{smp}.extendedFrags.fastq",
		R1="FilteredFastq/{smp}.notCombined_1.fastq",
		R2="FilteredFastq/{smp}.notCombined_2.fastq"
	log:
		"logs/logsFLASH/{smp}_flash.log"
	shell:
		"""
		scif --quiet run flash2 -M 250 {input.pairs_R1} {input.pairs_R2} -o {params.prefix} &> {log}
		"""
# Create a single file for merged and windows reads
rule Concatenate_Widows_And_Merged:
	input:
		Merged="FilteredFastq/{smp}.extendedFrags.fastq",
		WI="FilteredFastq/filtered_filter2_{smp}_widows.fastq"
	output:
		"FilteredFastq/filtered_{smp}_merged_widows.fastq"
	shell:
		"""
		cat {input.Merged} {input.WI} > {output}
		"""
#At the end of this first cleaning step, three fastq files are produced. A merger reads and widows reads file (... merged_widows.fastq) and two files for unconcanate read  (...notCombined_1.fastq, ...notCombined_2.fastq")

R1list=expand("FilteredFastq/{smp}.notCombined_1.fastq",smp=SAMPLES)
R2list=expand("FilteredFastq/{smp}.notCombined_2.fastq",smp=SAMPLES)
PElist=expand("FilteredFastq/filtered_{smp}_merged_widows.fastq",smp=SAMPLES)

#First step of de-novo"co-meta_assembly" with Megahit (Bruijn graph).
rule Megahit_Assembly:
	input:
		R1s = R1list,
		R2s = R2list,
		PEs = PElist
	params:
		prefix="Assembly_results",
		commaR1s = ",".join(R1list),
		commaR2s = ",".join(R2list),
		commaPEs = ",".join(PElist)
	output:
		"Assembly_results/{RUN}.contigs.fa"
	log:
		"logs/logsAssembly/Megahit_{RUN}.log"
	threads: 10
	shell:
		"""
		scif --quiet run megahit -t {threads}  -m 180e9  -1 {params.commaR1s} -2 {params.commaR2s} -r {params.commaPEs} -o {params.prefix} --out-prefix {RUN} --continue  2> {log}
		"""

#Second step of assembly (re-assembly) using overlaps between reads.
rule Cap3_Assembly:
	input:
		"Assembly_results/{RUN}.contigs.fa"
	output:
		ass="Assembly_results/{RUN}_Assembly_results",
		cont="Assembly_results/{RUN}.contigs.fa.cap.contigs",
		sig="Assembly_results/{RUN}.contigs.fa.cap.singlets",
	log:
		"logs/logsAssembly/CAP3_{RUN}.log"
	shell:
		"""
		scif --quiet run cap3 {input}>{output.ass} 2> {log}
		"""
#Production of a single file of contigs.
rule Merge_Mega_cap_contigs:
	input:
		cont="Assembly_results/{RUN}.contigs.fa.cap.contigs",
		sig="Assembly_results/{RUN}.contigs.fa.cap.singlets",
	output:
		"Assembly_results/{RUN}_contigs_assembly_results.fa"
	shell:
		"""
		cat {input.cont} {input.sig} > {output}
		"""


rule Assembly_informations:
	input:
		"Assembly_results/{RUN}_contigs_assembly_results.fa"
	output:
		"logs/logsAssembly/assembly_stats_{RUN}.txt"
	shell:
		"""
		scif --quiet run bioawk -c fastx '{{print $name, length($seq)}}' {input}  > {output}
		"""

#Indexing contigs file for later step of mapping.
rule Index_Assembly:
	input:
		"Assembly_results/{RUN}_contigs_assembly_results.fa"
	output:
		"Assembly_results/{RUN}_contigs_assembly_results.fa.bwt"
	shell:
		"""
		scif --quiet run bwa index {input}
		"""
#Mapping of all cleaned reads on the "Metagenome" (all contigs).
rule Map_On_Assembly:
	input:
		R1_unflash="FilteredFastq/{smp}.notCombined_1.fastq",
		R2_unflash="FilteredFastq/{smp}.notCombined_2.fastq",
		flash="FilteredFastq/filtered_{smp}_merged_widows.fastq",
		CONTIGS = "Assembly_results/{RUN}_contigs_assembly_results.fa",
		INDEX = "Assembly_results/{RUN}_contigs_assembly_results.fa.bwt",
	output:
		unflash="MappingOnAssembly/{smp}_unflash_on_{RUN}.bam",
		flash="MappingOnAssembly/{smp}_flash_on_{RUN}.bam",
		rmdupunflash="MappingOnAssembly/rmdup/{smp}_rmdupunflash_on_{RUN}.bam",
		rmdupflash="MappingOnAssembly/rmdup/{smp}_rmdupflash_on_{RUN}.bam",
		merged="MappingOnAssembly/{smp}_on_{RUN}.bam",
	threads: 5
	shell:
		"""
		scif --quiet run bwa mem -t {threads} {input.CONTIGS} {input.R1_unflash} {input.R2_unflash}| samtools view -b - > {output.unflash}
		scif --quiet run bwa mem -t {threads} {input.CONTIGS} {input.flash} | samtools view -b - > {output.flash}
		scif --quiet run samtools rmdup -S  {output.unflash} - > {output.rmdupunflash}
		scif --quiet runsamtools rmdup -S  {output.flash} - > {output.rmdupflash}
		scif --quiet runsamtools merge {output.merged} {output.rmdupflash} {output.rmdupunflash}
		"""

#Sort bam files from the mapping.
rule Sort_bam_mapped:
	input:
		"MappingOnAssembly/{smp}_on_{RUN}.bam",
	output:
		"MappingOnAssembly/sort_{smp}_on_{RUN}.bam"
	shell:
		"""
		scif --quiet run samtools sort {input} > {output}
		"""

#Use the mapping to determine the insert sizes
rule get_insert_size_metric:
	input:
		"MappingOnAssembly/sort_{smp}_on_{RUN}.bam"
	output:
		metrics="logs/insert_size/{smp}_insert_size_metrics_{RUN}.txt",
		histo="logs/insert_size/{smp}_insert_size_histo_{RUN}.pdf",
	shell:
		"""
		scif --quiet run picard CollectInsertSizeMetrics I={input}  O={output.metrics} H={output.histo}
		"""

#Extract some statistic information from the mapping.
rule Mapping_information:
	input:
		"MappingOnAssembly/{smp}_on_{RUN}.bam",
	output:
		"logs/{smp}_stats_mapping_assembly_{RUN}.txt",
	shell:
		"""
		scif --quiet run samtools flagstat {input} -o {output};
		"""

#Qantify the contigs coverage for each input samples.
rule Quantify_contigs_coverage:
	input:
		"MappingOnAssembly/{smp}_on_{RUN}.bam"
	output:
		mapped= "CountsMapping/{smp}_counts_contigs_{RUN}.mat",
		Unmapped="CountsMapping/{smp}_counts_unmapped_{RUN}.mat",

	shell:
		"""
		scif --quiet run samtools view -F 0x904 {input}| cut -f 3 | sort | uniq -c - > {output.mapped};
		scif --quiet run samtools view -f 0x4 {input}| cut -f 1 | sort | uniq -c - > {output.Unmapped};
		"""
#Extract some statistic information from contigs coverage quantification.
rule log_Quantify_contigs_coverage:
	input:
		mapped= "CountsMapping/{smp}_counts_contigs_{RUN}.mat",
		Unmapped="CountsMapping/{smp}_counts_unmapped_{RUN}.mat",
	output:
		"logs/logs_coverage/{smp}_coverage_{RUN}.txt"
	shell:
		"""
		awk '{{ sum+=$1 }} END {{print "Mapped:" sum }}' {input.mapped} >> {output}
		wc -l {input.mapped}  >> {output}
		awk '{{ sum+=$1 }}END {{print "UnMapped:" sum }}' {input.Unmapped} >> {output}
		wc -l {input.Unmapped}  >> {output}
		"""

# Extraction of reads which are not involved in the  contigs construction (not mapped to the metagenome).
rule Extract_filtered_Umapped_on_contigs:
	input:
		rmdupunflash="MappingOnAssembly/rmdup/{smp}_rmdupunflash_on_{RUN}.bam",
		rmdupflash="MappingOnAssembly/rmdup/{smp}_rmdupflash_on_{RUN}.bam",
	output:
		unflash="Unmapped/{smp}_unflash_unmapped_{RUN}.fastq",
		flash="Unmapped/{smp}_flash_unmapped_{RUN}.fastq",
	shell:
		"""
		samtools view -b -hf 0x4 {input.rmdupunflash} |scif --quiet run  samtools bam2fq - | scif --quiet run seqtk seq -A - > {output.unflash};
		samtools view -b -hf 0x4 {input.rmdupflash} | scif --quiet run  samtools bam2fq - | scif --quiet run seqtk seq -A - > {output.flash};
		"""


#First step of alignment of contigs  with diamond on viral protein database.
rule Blast_contigs_on_nr_vir:
	input:
		"Assembly_results/{RUN}_contigs_assembly_results.fa",
	output:
		"Blast_nr_vir_results/Contigs_{RUN}.blast_nr_vir_results.tsv",
	params:
		blastDBpath=base_nr_vir,
		basetaxoDBpath=base_taxo,
	threads: 5
	shell:
		"""
		scif --quiet run diamond blastx  -b 5 -c 1 -d {params.blastDBpath} --taxonlist 10239 --more-sensitive --query {input} --max-hsps 1 --max-target-seqs 1  --taxonmap {params.basetaxoDBpath} -f 6 qseqid sseqid qlen slen length qstart qend sstart send qcovhsp pident evalue bitscore staxids --out {output};
		"""

rule Remove_poor_qual_Blast_unmapped_on_nr_vir:
	input:
		"Blast_nr_vir_results/Contigs_{RUN}.blast_nr_vir_results.tsv",
	output:
		"Blast_nr_vir_results/Contigs_{RUN}.blast_nr_vir_results_cleaned.tsv",
	shell:
		"""
		awk '$3 >= 150 {{ print }}' {input} | awk '$5 >= 75 {{ print }}' > {output}
		"""

#Extraction of contigs with a hit on the viral database
rule Extract_hits_contigs_on_nr_vir:
	input:
		contigs="Assembly_results/{RUN}_contigs_assembly_results.fa",
		blast_contigs="Blast_nr_vir_results/Contigs_{RUN}.blast_nr_vir_results_cleaned.tsv",
	output:
		"Blast_nr_vir_hits_seq/hits_contigs_{RUN}.fa",
	params:
		scriptdir+"extract_blast_hits.py"
	shell:
		"""
		python {params} {input.blast_contigs} {input.contigs} {output}
		"""

#Second step of alignment of contigs with diamond on generalist protein database.
rule Blast_contigs_on_nr:
	input:
		"Blast_nr_vir_hits_seq/hits_contigs_{RUN}.fa"
	output:
		"Blast_nr_results/Contigs_{RUN}.blast_nr_results.tsv"
	params:
		blastDBpath=base_nr,
		basetaxoDBpath=base_taxo,
	shell:
		"""
		scif --quiet run diamond  blastx -d {params.blastDBpath} --sensitive --query {input} --max-hsps 1 --max-target-seqs 5  --taxonmap {params.basetaxoDBpath} -f 6 qseqid sseqid qlen slen length qstart qend sstart send qcovhsp pident evalue bitscore staxids --out {output};
		"""


# Extract some stats (logs) from the alignement steps
rule log_diamond:
	input:
		blastcontigs_nrvir="Blast_nr_vir_results/Contigs_{RUN}.blast_nr_vir_results.tsv",
		blastcontigs_nr="Blast_nr_results/Contigs_{RUN}.blast_nr_results.tsv"
	output:
		"logs/logs_diamond/{smp}_diamond.txt",
	shell:
		"""
		awk '{{print $1}}' {input.blastcontigs_nrvir} | sort -u | wc -l >> {output}
		awk '{{print $1}}' {input.blastcontigs_nr} | sort -u | wc -l >> {output}
		"""

unflash_list=expand("Blast_nr_results/{smp}_unflash_unmapped_"+RUN+".blast_nr_results.tsv", smp=SAMPLES)
flash_list=expand("Blast_nr_results/{smp}_flash_unmapped_"+RUN+".blast_nr_results.tsv", smp=SAMPLES)

#Production of a single file for all alignement results.
rule Join_seq_acc_taxo_nr :
	input:
		contigs="Blast_nr_results/Contigs_{RUN}.blast_nr_results.tsv",
	output:
		"Taxonomy/Seq_hits_info_{RUN}.csv",
	params:
	shell:
		"""
		cat {input.contigs} | awk -F'\t' '$14!=""' | sort -u -k1,1 | sed "s/;/\t/g" | awk '{{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$(NF)}}'| sed "s/ /\t/g" > {output}
		"""

#Use the taxonomic_id from the alignement to retreive the complete lineage.
rule get_nr_lineage_from_taxids:
	input:
		"Taxonomy/Seq_hits_info_{RUN}.csv"
	output:
		"Taxonomy/lineage_{RUN}.csv"
	params: scriptdir+"get_rank.py"
	shell:
		"""
		sort -u -k1,1  {input} |sed -e 's/;/\t/g'| awk '{{print $NF}}'| sort -u | sed '/^$/d' | paste -s -d,| python {params} {output}
		"""

count_contigs_list=expand("CountsMapping/{smp}_counts_contigs_"+RUN+".mat",smp=SAMPLES),
count_unmapped_list=expand("CountsMapping/{smp}_counts_unmapped_"+RUN+".mat",smp=SAMPLES),
count_list=count_contigs_list+count_unmapped_list

#Crosses, alignment, taxonomy, and contig coverage file to produce a single file of results
rule Build_array_coverage_nr:
	input:
		blast_info="Taxonomy/Seq_hits_info_{RUN}.csv",
		lineage="Taxonomy/lineage_{RUN}.csv",
		count_list=count_list
	output:
		lineage="Coverage/lineage_stats_{RUN}.csv",
		by_seq="Coverage/stats_by_seq_{RUN}.csv",
	params:
		script = scriptdir+"build_tables_coverage.py",
		countdir = "CountsMapping/",
	shell:
		"""
		python {params.script} {input.blast_info} {params.countdir} {input.lineage} {output.lineage} {output.by_seq}
		"""
#Use different database and custom file to complete taxonomic informations
rule complete_taxo:
	input:
		lineage="Coverage/lineage_stats_{RUN}.csv",
		by_seq="Coverage/stats_by_seq_{RUN}.csv",
	output:
		lineage="Coverage/abund_by_otu_{RUN}.csv",
		by_seq="Coverage/abund_by_contigs_{RUN}.csv",
	params:
		script = scriptdir+"complete_taxo.py",
		ref_taxo=ref_taxo,
		ref_taxo1=ref_taxo1,
		ICTV_taxo=ICTV_taxo,
		host_taxo=host_taxo,
	shell:
		"""
		python {params.script} {input.lineage} {input.by_seq} {params.ref_taxo} {params.ref_taxo1} {params.ICTV_taxo} {params.host_taxo} {output.lineage} {output.by_seq}
		"""

#Use logs from different rules to produce a summary off the analysis.
rule Create_logs_report:
	input:
		lineage="Coverage/abund_by_otu_{RUN}.csv",
		by_seq="Coverage/abund_by_contigs_{RUN}.csv",
	output:
		"logs/stats_run_{RUN}.csv"
	params:
		files=datadir,
		ext=ext_R1+ext,
		script=scriptdir+"create_results_doc.py"
	shell:
		"""
		python {params.script} {params.files} {params.ext} {input.lineage}  {input.by_seq} {output}
		"""

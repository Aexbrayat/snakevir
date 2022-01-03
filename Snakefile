import re
import sys
import os
from os.path import join
import shutil


def message(mes):
	sys.stderr.write("|---- " + mes + "\n")

def errormes(mes):
	sys.stderr.write("| ERROR ----" + mes + "\n")

configfile : "config.yaml"
datadir = config["fastq"]
scriptdir = config["Scripts"]
host_db = config["host_db"]
rRNA_bact=config["rRNA_bact"]
rRNA_host=config["rRNA_host"]
base_nr = config["base_nr"]
base_taxo = config["base_taxo"]
base_nt = config["base_nt"]
base_taxo_nt=config["base_taxo_nt"]
ext=config["ext"]
ext_R1=config["ext_R1"]
ext_R2=config["ext_R2"]
threads_default= config["threads_default"]
threads_Map_On_host= config["threads_Map_On_host"]
threads_Map_On_bacteria= config["threads_Map_On_bacteria"]
threads_Megahit_Assembly= config["threads_Megahit_Assembly"]
threads_Map_On_Assembly= config["threads_Map_On_Assembly"]
threads_Blast_contigs_on_nr= config["threads_Blast_contigs_on_nr"]
threads_Blast_contigs_on_nt= config["threads_Blast_contigs_on_nt"]

logclust  = os.getcwd()+"/logsclust/"

message(str(logclust))

READS, =  glob_wildcards(datadir+"{readfile}"+ext)
SAMPLES, = glob_wildcards(datadir+"{sample}"+ ext_R1+ext)


SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
snakemake.utils.makedirs("cluster_log/")


NBSAMPLES = len(SAMPLES)
NBREADS = len(READS)
RUN = config["run"]
#message(str(READS))
message(str(NBSAMPLES)+" samples  will be analysed")
message(str(len(READS))+" fastq files  will be processed")
message("Run name: "+RUN)
if NBREADS != 2*NBSAMPLES:
	errormes("Please provide two reads file per sample")
	sys.exit()

rule All:
	input:
		expand("cutadaptfiles/{readfile}.trimmed.fastq", readfile=READS),
		expand("logs/logs_contaminent/Stats_contaminent_{smp}.txt", smp=SAMPLES),
		"logs/logsAssembly/"+RUN+"_assembly_stats.txt",
		"Assembly_results/"+RUN+"_contigs_assembly_results.fa.bwt",
		expand("logs/insert_size/{smp}_insert_size_metrics_"+RUN+".txt", smp=SAMPLES),
		expand("logs/{smp}_"+RUN+"_stats_mapping_assembly.txt", smp=SAMPLES),
		expand("CountsMapping_raw/{smp}_raw_counts_contigs_"+RUN+".mat", smp=SAMPLES),
		expand("CountsMapping/{smp}_counts_contigs_"+RUN+".mat", smp=SAMPLES),
		expand("logs/logs_coverage/{smp}_coverage_"+RUN+".txt", smp=SAMPLES),
		"Taxonomy/lineage_"+RUN+".csv",
		"Coverage/count_table_contigs_"+RUN+".csv",
		"Coverage/count_contigs_raw_"+RUN+".csv",
		"Coverage/lineage_cor_"+RUN+".csv",
		expand("logs/logs_coverage_raw/{smp}_coverage_"+RUN+".txt", smp=SAMPLES),
		"Assembly_results/"+RUN+"_viral_contigs.fa",
		"Taxonomy_nt/Seq_hit_nt_"+RUN+".csv",
		"Blast_nt_results/Contigs_"+RUN+".blast_nt_results.tsv",
		"Taxonomy_nt/lineage_nt_"+RUN+".csv",
		"intergreted_vir_check_"+RUN+".csv",
		"Taxonomy/hosts_lineage_"+RUN+".csv",
		"logs/stats_run_"+RUN+".csv",
		dynamic(expand("depthvir/{num}/{smp}_on_"+RUN+".cov", smp=SAMPLES, num="{num}")),
		dynamic("depthvir/coverage_{num}.cov"),
		dynamic("depthvir/stats_coverage_{num}.cov"),
		dynamic("depthvir/range_10x_coverage_{num}.cov"),
		"results/range_10x_coverage.cov",

#Remove sequencing adapters based on 5' and 3' sequences (A3 and A5 in the config.yaml file)
#-n Trimming n adapter from each read ; -g Regular 5’ adapter; -a Regular 3’ adapter ; -j num core 0=auto number of available ; -o output
rule Remove_sequencing_adapters:
	input:
		datadir+"{readfile}"+ext ,
	output:
		"cutadaptfiles/{readfile}.trimmed.fastq",
	log:
		log="logs/logscutadapt/{readfile}_cut1.log",
	params:
		A3 = config["A3"],
		A5 = config["A5"],
	shell:
		"""
		cutadapt -n 10 -g {params.A5} -a {params.A3} -j 0 --overlap 15  {input} -o {output} &> {log.log}
		"""

# Quality trimming based on the --quality-cutoff (low-quality ends) -q parameter to change the cutoff
#and --minimum-length -m option to Discard processed reads that are shorter than the length specified.
# -q quality-cutoff 5',3' ; -j num core 0=auto number of available ; m minimum-length ; -o output
rule Quality_trimming:
	input:
		"cutadaptfiles/{readfile}.trimmed.fastq"
	output:
		"cutadaptfiles/{readfile}.clean.fastq"
	log:
		"logs/logscutadapt/{readfile}_cut2.log"
	benchmark:
		"benchmarks/{readfile}.cut2.benchmark.txt"
	shell:
		"""
		cutadapt -q 30,30 -j 0 -m 40 -o {output} {input} &> {log}
		"""

# Resynchronize 2 fastq or fastq.gz files (R1 and R2) after they have been trimmed and cleaned.
rule Repair_Pairs:
	input:
		R1="cutadaptfiles/{smp}"+ext_R1+".clean.fastq",
		R2="cutadaptfiles/{smp}"+ext_R2+".clean.fastq"
	output:
		R1="cutadaptfiles/{smp}"+ext_R1+".clean.fastq_pairs_R1.fastq",
		R2="cutadaptfiles/{smp}"+ext_R2+".clean.fastq_pairs_R2.fastq",
		WI="cutadaptfiles/{smp}"+ext_R1+".clean.fastq_singles.fastq"
	params:
		scriptdir+"fastqCombinePairedEnd.py"
	log:
		"logs/logsRepairspairs/{smp}_repair.log"
	benchmark:
		"benchmarks/{smp}.cut2.benchmark.txt"
	shell:
		"""
		python2 {params} {input.R1} {input.R2}  2> log
		"""

 # Host-homologous sequence cleaning by mapping (bwa) on ribosomal sequences.
rule Map_On_host:
	input:
		host = config["rRNA_host"],
		hostindex = config["rRNA_host"]+".bwt",
		R1="cutadaptfiles/{smp}"+ext_R1+".clean.fastq_pairs_R1.fastq",
		R2="cutadaptfiles/{smp}"+ext_R2+".clean.fastq_pairs_R2.fastq",
		WI="cutadaptfiles/{smp}"+ext_R1+".clean.fastq_singles.fastq"
	output:
		bam_pairs="HostMapping/{smp}_dipteria_pairs.bam",
		bam_WI="HostMapping/{smp}_dipteria_widows.bam",
		sam_pairs="HostMapping/unmapped_{smp}_dipteria_pairs.sam",
		sam_WI="HostMapping/unmapped_{smp}_dipteria_widows.sam"
	benchmark:
		"benchmarks/{smp}.Map_Pairs_On_host.benchmark.txt"
	log:
		"logs/logsMapHost/{smp}_bwa_pairs_on_dipt.log"
	threads: threads_Map_On_host
	shell:
		"""
		if [ -s {input.R1} ] && [ -s {input.R2} ]
		then
			bwa mem -t {threads} {input.host} {input.R1} {input.R2} 2> {log} |samtools view -b - | tee {output.bam_pairs} | samtools view -f 0x4 -> {output.sam_pairs}
		else
			echo "{input.R1} or {input.R2} is empty."
			touch {output.bam_pairs}
			touch {output.sam_pairs}
		fi

		if [ -s {input.WI} ]
		then
		bwa mem -t {threads} {input.host} {input.WI} 2> {log} | samtools view -b - | tee {output.bam_WI} | samtools view -f 0x4 -> {output.sam_WI}
		else
			echo "{input.WI} is empty."
			touch {output.bam_WI}
			touch {output.sam_WI}
		fi
		"""

# Extact samples squences with with no homologous sequences with the host.
rule Extract_Unmapped_host_Reads:
	input:
		sam_pairs="HostMapping/unmapped_{smp}_dipteria_pairs.sam",
		sam_WI="HostMapping/unmapped_{smp}_dipteria_widows.sam"
	output:
		pairs_R1="FilteredFastq/filtered_diptera_{smp}_R1.fastq",
		pairs_R2="FilteredFastq/filtered_diptera_{smp}_R2.fastq",
		WI="FilteredFastq/filtered_diptera_{smp}_widows.fastq"
	benchmark:
		"benchmarks/{smp}.Extract_Unmapped_diptera_Reads.benchmark.txt"
	shell:
		"""
		if [ -s {input.sam_pairs} ]
		then
			picard SamToFastq VALIDATION_STRINGENCY=SILENT I={input.sam_pairs} F={output.pairs_R1} F2={output.pairs_R2}
		else
			echo "{input.sam_pairs} or {input.sam_WI} is empty."
			touch {output.pairs_R1}
			touch {output.pairs_R2}
		fi
		if [ -s {input.sam_WI} ]
		then
			picard SamToFastq VALIDATION_STRINGENCY=SILENT I={input.sam_WI} F={output.WI}
		else
			echo "{input.sam_WI} is empty."
			touch {output.WI}
		fi
		"""

# Bacterial-homologous sequences cleaning by mapping (bwa) on ribosomal sequences.
rule Map_On_bacteria:
	input:
		host = config["rRNA_bact"],
		hostindex = config["rRNA_bact"]+".bwt",
		pairs_R1="FilteredFastq/filtered_diptera_{smp}_R1.fastq",
		pairs_R2="FilteredFastq/filtered_diptera_{smp}_R2.fastq",
		WI="FilteredFastq/filtered_diptera_{smp}_widows.fastq"
	output:
		bam_pairs="HostMapping/{smp}_bacteria_pairs.bam",
		bam_WI="HostMapping/{smp}_bacteria_widows.bam",
		sam_pairs="HostMapping/unmapped_{smp}_bacteria_pairs.sam",
		sam_WI="HostMapping/unmapped_{smp}_bacteria_widows.sam"
	benchmark:
		"benchmarks/{smp}.Map_Pairs_On_bact.benchmark.txt"
	log:
		"logs/logsMapBact/{smp}_bwa_pairs_on_bact.log"
	threads: threads_Map_On_bacteria
	shell:
		"""

		if [ -s {input.pairs_R1} ] && [ -s {input.pairs_R2} ]
		then
			bwa mem -t {threads} {input.host} {input.pairs_R1} {input.pairs_R2} 2> {log} | samtools view -b -  | tee {output.bam_pairs} |  samtools view -f 0x4 -> {output.sam_pairs}
		else
			echo "{input.pairs_R1} or {input.pairs_R2} is empty."
			touch {output.bam_pairs}
			touch {output.sam_pairs}
		fi

		if [ -s {input.WI} ]
		then
		bwa mem -t {threads} {input.host} {input.WI} 2> {log} | samtools view -b - | tee {output.bam_WI} | samtools view -f 0x4 -> {output.sam_WI}
		else
			echo "{input.WI} is empty."
			touch {output.bam_WI}
			touch {output.sam_WI}
		fi
		"""

# Extact samples squences with with no homologous sequences with bacteria.
rule Extract_Unmapped_bact_Reads:
	input:
		sam_pairs="HostMapping/unmapped_{smp}_bacteria_pairs.sam",
		sam_WI="HostMapping/unmapped_{smp}_bacteria_widows.sam"
	output:
		pairs_R1="FilteredFastq/filtered_bacteria_{smp}_R1.fastq",
		pairs_R2="FilteredFastq/filtered_bacteria_{smp}_R2.fastq",
		WI="FilteredFastq/filtered_bacteria_{smp}_widows.fastq"
	benchmark:
		"benchmarks/{smp}.Extract_Unmapped_bacteria_Reads.benchmark.txt"
	shell:
		"""
		if [ -s {input.sam_pairs} ]
		then
			picard  SamToFastq VALIDATION_STRINGENCY=SILENT I={input.sam_pairs} F={output.pairs_R1} F2={output.pairs_R2}
		else
			echo "{input.sam_pairs} is empty."
			touch {output.pairs_R1}
			touch {output.pairs_R2}
		fi
		if [ -s {input.sam_WI} ]
		then
			picard SamToFastq  VALIDATION_STRINGENCY=SILENT I={input.sam_WI} F={output.WI}
		else
			echo "{input.sam_WI} is empty."
			touch {output.WI}
		fi
		"""
 # Extract some stats (logs) from the cleaning steps
rule log_map_conta:
	input:
		host_pairs_R1="FilteredFastq/filtered_diptera_{smp}_R1.fastq",
		host_pairs_R2="FilteredFastq/filtered_diptera_{smp}_R2.fastq",
		host_WI="FilteredFastq/filtered_diptera_{smp}_widows.fastq",
		bact_pairs_R1="FilteredFastq/filtered_bacteria_{smp}_R1.fastq",
		bact_pairs_R2="FilteredFastq/filtered_bacteria_{smp}_R2.fastq",
		bact_WI="FilteredFastq/filtered_bacteria_{smp}_widows.fastq"
	output:
		"logs/logs_contaminent/Stats_contaminent_{smp}.txt",
	shell:
		"""
		awk '{{s++}}END{{print "Host_pair_R1 : " s/4}}' {input.host_pairs_R1} >> {output};
		awk '{{s++}}END{{print "Host_pairs_R2 : " s/4}}' {input.host_pairs_R2} >> {output};
		awk '{{s++}}END{{print "Host_pair_wi : " s/4}}' {input.host_WI} >> {output};
		awk '{{s++}}END{{print "bact_pair_R1 : " s/4}}' {input.bact_pairs_R1} >> {output};
		awk '{{s++}}END{{print "bact_pairs_R2 : " s/4}}' {input.bact_pairs_R2} >> {output};
		awk '{{s++}}END{{print "bact_pair_wi : " s/4}}' {input.bact_WI} >> {output};
		"""

# Merge paired-end reads , 3 files where created on with de merde pairs and two for the unmerged R1 and R2
rule Merge_Pairs_With_Flash:
	input:
		pairs_R1="FilteredFastq/filtered_bacteria_{smp}_R1.fastq",
		pairs_R2="FilteredFastq/filtered_bacteria_{smp}_R2.fastq",
	params:
		prefix="FilteredFastq/{smp}",
		flash=scriptdir+"flash"
	output:
		ext="FilteredFastq/{smp}.extendedFrags.fastq",
		R1="FilteredFastq/{smp}.notCombined_1.fastq",
		R2="FilteredFastq/{smp}.notCombined_2.fastq"
	benchmark:
		"benchmarks/{smp}.Merge_Pairs_With_Flash.benchmark.txt"
	log:
		"logs/logsFLASH/{smp}_flash.log"
	shell:
		"""
		flash -M 250 {input.pairs_R1} {input.pairs_R2} -o {params.prefix} &> {log}
		"""

# Create a single file for merge and windows reads
rule Concatenate_Widows_And_Merged:
	input:
		Merged="FilteredFastq/{smp}.extendedFrags.fastq",
		WI="FilteredFastq/filtered_bacteria_{smp}_widows.fastq"
	output:
		"FilteredFastq/filtered_{smp}_merged_widows.fastq"
	benchmark:
		"benchmarks/{smp}.Concatenate_Widows_And_Merged.benchmark.txt"
	shell:
		"""
		cat {input.Merged} {input.WI} > {output}
		"""


#At the end of this first cleaning step, three fastq files are produced. A merger reads and widows reads file (... merged_widows.fastq) and two files for unconcanate read  (...notCombined_1.fastq, ...notCombined_2.fastq")

R1list=expand("FilteredFastq/{smp}.notCombined_1.fastq",smp=SAMPLES)
R2list=expand("FilteredFastq/{smp}.notCombined_2.fastq",smp=SAMPLES)
PElist=expand("FilteredFastq/filtered_{smp}_merged_widows.fastq",smp=SAMPLES)

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
	threads: threads_Megahit_Assembly
	shell:
		"""
		touch Assembly_results/{wildcards.RUN}.contigs.fa
		megahit -t {threads}  --mem-flag 2 -m 1 -1 {params.commaR1s} -2 {params.commaR2s} -r {params.commaPEs} -o Assembly_results_tmp --out-prefix {wildcards.RUN} --continue  2> {log}
		rm -r Assembly_results
		mv Assembly_results_tmp Assembly_results
		"""

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
		cap3 {input}>{output.ass} 2> {log}
		"""

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
		"logs/logsAssembly/{RUN}_assembly_stats.txt",
	shell:
		"""
		bioawk -c fastx '{{ print $name, length($seq) }}' {input} | sed 's/      / /g' > {output}
		"""

rule Index_Assembly:
	input:
		"Assembly_results/{RUN}_contigs_assembly_results.fa"
	output:
		"Assembly_results/{RUN}_contigs_assembly_results.fa.bwt"
	shell:
		"""
		bwa index {input} > {output}
		"""

rule Map_On_Assembly:
	input:
		pairs_R1="FilteredFastq/filtered_bacteria_{smp}_R1.fastq",
		pairs_R2="FilteredFastq/filtered_bacteria_{smp}_R2.fastq",
		WI="FilteredFastq/filtered_bacteria_{smp}_widows.fastq",
		CONTIGS = "Assembly_results/{RUN}_contigs_assembly_results.fa",
		INDEX = "Assembly_results/{RUN}_contigs_assembly_results.fa.bwt",
	output:
		raw="MappingOnAssembly/raw_{smp}_on_{RUN}.bam",
		raw_wi="MappingOnAssembly/raw_wi_{smp}_on_{RUN}.bam",
		rmdup="MappingOnAssembly/rmdup_{smp}_on_{RUN}.bam",
		rmdup_wi="MappingOnAssembly/rmdup_wi_{smp}_on_{RUN}.bam",
		merged="MappingOnAssembly/{smp}_on_{RUN}.bam",
		mergedraw="MappingOnAssembly/{smp}_raw_on_{RUN}.bam",
		logs_pair="logs/logsDuplicates/{smp}_duplicates_pairs_{RUN}.txt",
		logs_wi="logs/logsDuplicates/{smp}_duplicates_wi_{RUN}.txt",
		logs="logs/{smp}_{RUN}_stats_mapping_assembly.txt",
	threads: threads_Map_On_Assembly
	shell:
		"""
		if [ -s {input.pairs_R1} ]
		then
			bwa mem -t {threads} {input.CONTIGS} {input.pairs_R1} {input.pairs_R2}| samtools fixmate -@ {threads} -m - -| samtools sort -@ {threads} -O BAM -|  tee {output.raw} | samtools markdup -@ {threads} -r - -f {output.logs_pair} {output.rmdup}
		else
			touch {output.raw} {output.rmdup}
		fi
		if [ -s {input.WI} ]
		then
			bwa mem -t {threads} {input.CONTIGS} {input.WI}| samtools fixmate -@ {threads} -m - -| samtools sort -@ {threads} -O BAM |  tee {output.raw_wi} | samtools markdup -@ {threads} -r - -f {output.logs_wi}  {output.rmdup_wi}
		else
			touch {output.raw_wi} {output.rmdup_wi}
		fi
		samtools merge -@ {threads} {output.merged} {output.rmdup} {output.rmdup_wi}
		samtools merge -@ {threads} {output.mergedraw} {output.raw} {output.raw_wi}
		samtools flagstat {output.merged} > {output.logs}
		"""

rule get_insert_size_metric:
	input:
		"MappingOnAssembly/raw_{smp}_on_{RUN}.bam",
	output:
		metrics="logs/insert_size/{smp}_insert_size_metrics_{RUN}.txt",
		histo="logs/insert_size/{smp}_insert_size_histo_{RUN}.pdf",
	shell:
		"""
		picard  CollectInsertSizeMetrics -M 0.0 -I {input}  -O {output.metrics} -H {output.histo}
		"""



rule Quantify_contigs_coverage_raw:
	input:
		"MappingOnAssembly/{smp}_raw_on_{RUN}.bam",
	output:
		mapped= "CountsMapping_raw/{smp}_raw_counts_contigs_{RUN}.mat",
		Unmapped="CountsMapping_raw/{smp}_raw_counts_unmapped_{RUN}.mat",
	shell:
		"""
		echo "qseqid "$(basename {input} _raw_on_{RUN}.bam)> {output.mapped}
		samtools view -F 2308 {input}| cut -f 3 | sort | uniq -c - | awk -F' ' '{{print $NF,$(NF-1)}}' >> {output.mapped};
		echo "qseqid "$(basename {input} _raw_on_{RUN}.bam)> {output.Unmapped}
		samtools view -f 0x4 {input}| cut -f 1 | sort | uniq -c - | awk -F' ' '{{print $NF,$(NF-1)}}' >> {output.Unmapped};
		"""


rule Quantify_contigs_coverage_rmdup:
	input:
		"MappingOnAssembly/{smp}_on_{RUN}.bam"
	output:
		mapped= "CountsMapping/{smp}_counts_contigs_{RUN}.mat",
		Unmapped="CountsMapping/{smp}_counts_unmapped_{RUN}.mat",
	shell:
		"""
		echo "qseqid "$(basename {input} _on_{RUN}.bam)> {output.mapped}
		samtools view -F 2308 {input}| cut -f 3 | sort | uniq -c - | awk -F' ' '{{print $NF,$(NF-1)}}' >> {output.mapped};
		echo "qseqid "$(basename {input} _on_{RUN}.bam)> {output.Unmapped}
		samtools view -f 0x4 {input}| cut -f 1 | sort | uniq -c - | awk -F' ' '{{print $NF,$(NF-1)}}' >> {output.Unmapped};
		"""


rule log_Quantify_contigs_coverage:
	input:
		mapped= "CountsMapping/{smp}_counts_contigs_{RUN}.mat",
		Unmapped="CountsMapping/{smp}_counts_unmapped_{RUN}.mat",
	output:
		"logs/logs_coverage/{smp}_coverage_{RUN}.txt"
	shell:
		"""
		awk '{{ sum+=$2 }} END {{print "Mapped:" sum }}' {input.mapped} >> {output}
		wc -l {input.mapped}  >> {output}
		awk '{{ sum+=$2 }}END {{print "UnMapped:" sum }}' {input.Unmapped} >> {output}
		wc -l {input.Unmapped}  >> {output}
		"""


rule Blast_contigs_on_nr:
	input:
		"Assembly_results/{RUN}_contigs_assembly_results.fa",
	output:
		"Blast_nr_results/Contigs_{RUN}.blast_nr_results_raw.tsv"
	params:
		blastDBpath=base_nr,
		basetaxoDBpath=base_taxo,
	threads: threads_Blast_contigs_on_nr
	shell:
		"""
		diamond blastx -b 10.0 -c 1 -p {threads}  -d {params.blastDBpath} --more-sensitive --query {input} --max-hsps 1 --max-target-seqs 5  --taxonmap {params.basetaxoDBpath} -f 6 qseqid sseqid qlen slen length qstart qend sstart send qcovhsp pident evalue bitscore staxids --out {output};
		"""

rule Remove_poor_qual_Blast_unmapped_on_nr_vir:
	input:
		"Blast_nr_results/Contigs_{RUN}.blast_nr_results_raw.tsv"
	output:
		"Blast_nr_results/Contigs_{RUN}.blast_nr_results.tsv",
	threads: 5
	shell:
		"""
		awk '$3 >= 150 {{ print }}' {input} | awk '$5 >= 75 {{ print }}' > {output}
		"""

rule log_diamond:
	input:
		blastcontigs_nrvir="Blast_nr_vir_results/Contigs_{RUN}.blast_nr_vir_results.tsv",
		blastcontigs_nr="Blast_nr_results/Contigs_{RUN}.blast_nr_results.tsv"
	output:
		"logs/logs_diamond/{smp}_{RUN}_diamond.txt",
	shell:
		"""
		awk '{{print $1}}' {input.blastcontigs_nrvir} | sort -u | wc -l >> {output}
		awk '{{print $1}}' {input.blastcontigs_nr} | sort -u | wc -l >> {output}
		"""
#
rule Join_seq_acc_taxo_nr :
	input:
		contigs="Blast_nr_results/Contigs_{RUN}.blast_nr_results.tsv",
	output:
		"Taxonomy/Seq_hits_info_{RUN}.csv",
	shell:
		"""
		cat {input.contigs} | awk -F'\t' '$14!=""' | sort -u -k1,1 | sed "s/;/\t/g" | awk '{{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$(NF)}}'| sed "s/ /\t/g" |sed -e 's/\t/,/g'| sed -e '1i\qseqid,sseqid,qlen,slen,length,qstart,qend,sstart,send,qcovhsp,pident,evalue,bitscore,tax_id' > {output}
		"""

rule get_nr_lineage_from_taxids:
	input:
		seq="Taxonomy/Seq_hits_info_{RUN}.csv",
		blst="Blast_nr_results/Contigs_{RUN}.blast_nr_results.tsv"
	output:
		lin="Taxonomy/lineage_{RUN}.csv",
		lin5="Taxonomy/lineage_5_hit_{RUN}.csv",
		lin5c="Taxonomy/lineage_5_hit_{RUN}_cor.csv",
		blst="Taxonomy/Blast_5_hit_{RUN}.csv",
		sort_seq_ids="tmp/sort_seq_ids_{RUN}.csv",

	params:
		getr=scriptdir+"get_rank.py",
		com=scriptdir+"complet_taxo_dic_v3.py",
		pickle_shi=scriptdir+"correc_taxo.pickle",
		pickle_cust=scriptdir+"custom_taxo.pickle"
	shell:
		"""
		sort -u -k1,1  {input.seq} |sed -e 's/;/,/g'| awk -F',' '{{print $NF}}'| sort -u | sed '/^$/d' | paste -s -d,| python {params.getr} {output.lin}
		cat {input.blst} | awk  '$14!=""' |sed "s/;/\t/g" | awk '{{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$(NF)}}'| sed "s/ /,/g" |sed -e 's/\t/,/g' |awk -F',' '{{print $NF}}'| sort -u | sed '/^$/d' | paste -s -d,  | python {params.getr} {output.lin5}
		cat {input.blst} | awk  '$14!=""' |sed "s/;/\t/g" | awk '{{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$(NF)}}'| sed "s/ /,/g" |sed -e 's/\t/,/g'| sed -e '1i\qseqid,sseqid,qlen,slen,length,qstart,qend,sstart,send,qcovhsp,pident,evalue,bitscore,tax_id' > {output.blst}
		python {params.com} {params.pickle_shi} {params.pickle_cust} {output.lin5} {output.lin5c}
		awk -F "," '{{print $NF, $1}}' {input.seq} | sort -n -k1 > {output.sort_seq_ids}

		"""




##count_contigs_list=expand("CountsMapping/{smp}_counts_contigs_"+RUN+".mat",smp=SAMPLES),
##count_unmapped_list=expand("CountsMapping/{smp}_counts_unmapped_"+RUN+".mat",smp=SAMPLES),
##count_list=count_contigs_list+count_unmapped_list

rule Build_array_coverage_nr:
	input:
		lin="Taxonomy/lineage_{RUN}.csv",
	output:
		by_seq="Coverage/count_table_contigs_{RUN}.csv",
	params:
		script = scriptdir+"build_count_table.py",
		countdir = "CountsMapping/",
	shell:
		"""
		python {params.script} {params.countdir} {output.by_seq}
		"""

rule Count_coverage_raw_nr:
	input:
		"Coverage/count_table_contigs_{RUN}.csv"
	output:
		by_seq="Coverage/count_contigs_raw_{RUN}.csv",
	params:
		script = scriptdir+"build_count_table.py",
		countdir = "CountsMapping_raw/",
	shell:
		"""
		python {params.script} {params.countdir} {output.by_seq}
		"""

rule complete_taxo:
	input:
		blast="Taxonomy/Blast_5_hit_{RUN}.csv",
		lin="Taxonomy/lineage_{RUN}.csv",
		lin5c="Taxonomy/lineage_5_hit_{RUN}_cor.csv",

	output:
		lineage1="Coverage/lineage_{RUN}_tmp.csv",
		lineage="Coverage/lineage_cor_{RUN}.csv",

	params:
		script_mult =scriptdir+"multihit.py",
		script_compt =scriptdir+"complet_taxo_dic_v3.py",
		pickle_shi=scriptdir+"correc_taxo.pickle",
		pickle_cust=scriptdir+"custom_taxo.pickle"
	shell:
		"""
		python {params.script_mult} {input.blast} {input.lin5c} {input.lin} tmp/lin.tmp
		cat tmp/lin.tmp |sed -e 's/,/;/g'|sed -e 's/\\t/,/g' > {output.lineage1}
		python {params.script_compt} {params.pickle_shi} {params.pickle_cust} {output.lineage1} {output.lineage}
		"""

rule extract_viral_ids:
	input:
		lineage="Coverage/lineage_cor_{RUN}.csv",
		sort_seq_ids="tmp/sort_seq_ids_{RUN}.csv",

	output:
		"Taxonomy/viral_contigs_ids_{RUN}.txt"
	shell:
		"""
		for i in `awk -F',' '{{print $1}}' {input.lineage}|  sort -n ` ; do grep $i {input.sort_seq_ids} | awk '{{print $2}}'> {output} ; done
		"""


rule extract_viral_contigs:
	input:
		contig = "Assembly_results/{RUN}_contigs_assembly_results.fa",
		viral_cont_IDS="Taxonomy/viral_contigs_ids_{RUN}.txt"
	output:
		viral_cont="Assembly_results/{RUN}_viral_contigs.fa",
	shell:
		"""
		perl -ne 'if(/^>(\S+)/){{$c=$i{{$1}}}}$c?print:chomp;$i{{$_}}=1 if @ARGV' {input.viral_cont_IDS} {input.contig} > {output.viral_cont};
		"""



rule log_Quantify_contigs_coverage_raw:
	input:
		mapped= "CountsMapping_raw/{smp}_raw_counts_contigs_{RUN}.mat",
		Unmapped="CountsMapping_raw/{smp}_raw_counts_unmapped_{RUN}.mat",
		viral_cont="Taxonomy/viral_contigs_ids_{RUN}.txt"
	output:
		"logs/logs_coverage_raw/{smp}_coverage_{RUN}.txt"
	shell:
		"""
		xargs -I @ grep -w -m 1 @ {input.mapped} | < {input.viral_cont}| awk '{{ sum+=$2 }} END {{print "Mapped:" sum }}' - >> {output}
		xargs -I @ grep -w -m 1 @ {input.mapped} | < {input.viral_cont} | wc -l >> {output}
		xargs -I @ grep -w -m 1 @ {input.Unmapped} | < {input.viral_cont}| awk '{{ sum+=$2 }} END {{print "Mapped:" sum }}' - >> {output}
		xargs -I @ grep -w -m 1 @ {input.Unmapped} | < {input.viral_cont} | wc -l >> {output}
		"""


rule Blast_contigs_on_nt:
	input:
		contig = "Assembly_results/{RUN}_contigs_assembly_results.fa",
		viral_cont="Assembly_results/{RUN}_viral_contigs.fa"
	output:
		blast="Blast_nt_results/Contigs_{RUN}.blast_nt_results.tsv",
	params:
		blastDBpath = base_nt
	threads: threads_Blast_contigs_on_nt
	shell:
		"""
		blastn -task blastn -db {params.blastDBpath} -query {input.viral_cont} -num_threads {threads}  -evalue 0.001  -max_hsps 1 -max_target_seqs 10 -outfmt "6 qseqid sseqid qlen slen length qstart qend sstart send qcovhsp pident evalue bitscore"  -out {output.blast}
		"""

rule extract_seq_acc:
	input:"Blast_nt_results/Contigs_{RUN}.blast_nt_results.tsv",
	output:"tmp/Contigs_{RUN}.blast_nt_results.tsv",
	shell:
		"""
		awk -F'|' '{{print $1,$(NF-1)}}' {input} | awk '{{$2="";print}}'| awk  '{{$2=substr($2,1, length($2)-2);print}}'| sort -u -k1,1 | sed 's/ /\t/g' |tee {output} | awk '{{print $2}}' >  "tmp/sseq_ids.txt";
		"""

rule extract_tax_ids:
	input:
		acc_ids= "tmp/Contigs_{RUN}.blast_nt_results.tsv",
	output:
		"tmp/tax_ids_{RUN}.tsv"
	params:
		basetaxoDBpath=base_taxo_nt,
	shell:
		"""
		for i in `cat {input.acc_ids}` ; do LC_ALL=C  look $i {params.basetaxoDBpath} | awk '{{print $1,$3}}' >> {output}; done || true
		"""

rule join_blst_tax:
	input:
		tax_ids="tmp/tax_ids_{RUN}.tsv",
		blst="Blast_nt_results/Contigs_{RUN}.blast_nt_results.tsv",
		blst_ids="tmp/Contigs_{RUN}.blast_nt_results.tsv",
	output:
		blst="Taxonomy_nt/Seq_hit_nt_{RUN}.csv"
	shell:
		"""
		awk 'NR==FNR{{a[$1]=$2;next}} ($2) in a{{print $0" "a[$2]}}'  {input.tax_ids}  {input.blst_ids} |awk '{{print $1 , $3}}'| sed 's/ /\t/g' > tmp/sseq_ids_tax_id.tsv
		cat {input.blst} | sort -u -k1,1| sed 's/ /\t/g' > tmp/blast.tsv
		awk 'NR==FNR{{a[$1]=$2;next}} ($1) in a{{print $0" "a[$1]}}' tmp/sseq_ids_tax_id.tsv  tmp/blast.tsv| sed 's/ /\t/g' | sed 's/\t/,/g' > {output.blst}
		"""
rule get_nt_lineage_from_taxids:
	input:
		"Taxonomy_nt/Seq_hit_nt_{RUN}.csv"
	output:
		lin="Taxonomy_nt/lineage_nt_{RUN}.csv",
		stat="Coverage_nt/stat_by_seq_nt_{RUN}.csv",
	params:
		scriptdir+"get_rank.py",
	shell:
		"""
		sort -u -k1,1  {input} |sed -e 's/;/,/g'| awk -F',' '{{print $NF}}'| sort -u | sed '/^$/d' | paste -s -d,| python {params} {output.lin}
		awk -F ',' 'NR==FNR{{a[$1]=$0;next}} ($NF) in a{{print $0","a[$NF]}}' {output.lin} {input} | sed -e '1i\qseqid,sseqid,qlen,slen,length,qstart,qend,sstart,send,qcovhsp,pident,evalue,bitscore,tax_id' > {output.stat}
		"""


rule intergreted_vir_check:
	input:
		stat_nt="Coverage_nt/stat_by_seq_nt_{RUN}.csv",
		stat="Taxonomy_nt/Seq_hit_nt_{RUN}.csv",
		lineage="Coverage/lineage_cor_{RUN}.csv",
	output:
		tmp="tmp/interg_virus_{RUN}.csv",
	shell:
		"""
		awk -F "," ' {{if ($17!="Viruses" && $17!="NA" && $10>25) print $0}}' {input.stat_nt} |sed -e 's/,/\t/g' > {output.tmp}
		"""

rule merge_vir_check:
	input:
		tmp="tmp/interg_virus_{RUN}.csv",
		stat="Taxonomy_nt/Seq_hit_nt_{RUN}.csv",
		lineage="Coverage/lineage_cor_{RUN}.csv",
	output:
		vir_ch="intergreted_vir_check_{RUN}.csv",

	shell:
		"""
		awk -F',' '{{print $1}}' {input.lineage}| grep -wf - {input.stat} | awk -F"[\t, ]" 'FNR==NR&&NR!=1{{a[$1]="YES"}}  FNR!=NR{{OFS=",";if(FNR!=1){{$(NF+1)=a[$1]}};print}}'  {input.tmp} - |awk -F',' '!$NF {{OFS=","; $NF="NO" }}1' >  {output.vir_ch}
		"""


rule split_viral_tax:
	input:
		"Coverage/lineage_cor_{RUN}.csv",
	output:
		dynamic("tmp/{RUN}_id_for_tax_{num}"),
	shell:
		"""
		 split -d -l 1000 {input}  {output}
		"""


rule get_host:
	input:
		acc_ids= "tmp/{RUN}_id_for_tax_{num}",
		by_seq="Taxonomy/Seq_hits_info_{RUN}.csv"
	output:
		"tmp/{RUN}_host_tax_{num}"
	params:
		script=scriptdir+"get_host.py",
		host_db=scriptdir+"virushostdb1_21.csv"
	shell:
		"""
		python {params.script} {input.acc_ids} {input.by_seq} {params.host_db} {output}
		"""

rule join_host_tax:
	input:
		dynamic("tmp/{RUN}_host_tax_{num}"),
	output:
		"Taxonomy/hosts_lineage_{RUN}.csv"
	shell:
		"""
		cat {input}  > {output}
		"""

rule Create_logs_report:
	input:
		lin="Taxonomy_nt/lineage_nt_{RUN}.csv",
		lineage="Coverage/lineage_cor_{RUN}.csv",
		by_seq="Taxonomy/Seq_hits_info_{RUN}.csv",
		count_t="Coverage/count_table_contigs_{RUN}.csv",
	output:
		"logs/stats_run_{RUN}.csv"
	params:
		files=datadir,
		ext=ext_R1+ext,
		script=scriptdir+"create_results_doc_new.py"
	shell:
		"""
		python {params.script} {params.files} {ext_R1} {ext_R2} {ext} {RUN} {input.by_seq}   {input.lineage}    {input.count_t}  {output}
		"""


rule split_depth:
	input:
		"MappingOnAssembly/{smp}_on_{RUN}.bam",
	output:
		"depth/{smp}_on_{RUN}.cov",
	shell:
		"""
		samtools depth -a {input} {output}
		"""


rule split_contigs_ids:
	input:
		"Taxonomy/viral_contigs_ids_{RUN}.txt",
	output:
		dynamic("tmp_spit_contigs/{RUN}_viral_contigs_ids_{num}"),
	shell:
		"""
		split -d -l 10000 {input} {output}
		"""

checkpoint Depth_vir:
	input:
		viral_cont_IDS="tmp_spit_contigs/{RUN}_viral_contigs_ids_{num}",
		depth="depth/{smp}_on_{RUN}.cov",
	output:
		"depthvir/{num}/{smp}_on_{RUN}.cov"
	shell:
		"""
		mkdir -p "depthvir/{wildcards.num}/"
		for i in `cat {input.viral_cont_IDS}` ; do rg -w $i {input.depth}  ; done > {output} || true
		"""

rule merge_depth:
	input:
	output:
		"depthvir/coverage_{num}.cov"
	params:
		dir=dynamic("depthvir/{num}/"),
		script=scriptdir+"merge_cov.py"
	shell:
		"""
		python {params.script} {params.dir} {output}
		"""

checkpoint stats_depth:
	input:
		"depthvir/coverage_{num}.cov"
	output:
		"depthvir/stats_coverage_{num}.cov"
	params:
		script=scriptdir+"stats_coverage.py"
	shell:
		"""
		python {params.script} {input} {output}
		"""

checkpoint range_depth:
	input:
		"depthvir/coverage_{num}.cov"
	output:
		"depthvir/range_10x_coverage_{num}.cov"
	params:
		script=scriptdir+"count_range_cov.py"
	shell:
		"""
		python {params.script} {input} {output}
		"""

rule results_depth:
	input:
		depthstat=dynamic("depthvir/stats_coverage_{num}.cov"),
		depth10x=dynamic("depthvir/range_10x_coverage_{num}.cov")
	output:
		stats_deph="results/stats_coverage.cov",
		range_10x="results/range_10x_coverage.cov"
	shell:
		"""
		echo "contig,min_coverage,max_coverage, mean_coverage,median_coverage" > {output.stats_deph}
		cat depthvir/stats_coverage_* >> {output.stats_deph}
		echo "contig,start,end" > {output.range_10x}
		cat  depthvir/range_10x_coverage_* >> {output.range_10x}
		"""

# vim: set ft=python:

# RNAseq Single-end workflow v1.0
# Hernan Lorenzi
# hernan.lorenzi@nih.gov
# Workflow requires to have a STAR index already available within the data/00ref directory
# Also, it is necessary to configure the config.yml file accordingly to include all metadatata required.

import os
import glob

configfile: "config/config.yaml"
samples = config["samples"].keys()
genome = config["reference"]["genome_file"]
annotation = config["reference"]["ensembl_gtf"]
starOverhang = config["star_db"]["sjdbOverhang"]

# Build directory structure
paths = ["data/00ref", "data/00RawData/","data/00reads", "data/00adapters",
        "results/01trim", "results/02abundant","results/03map_reads",
        "results/04dedup","results/05correlation", "results/05bigwig",
        "results/05counts","results/06fastqc_raw",
        "results/06fastqc_trim","results/07multiqc"
        ]
for path in paths:
    os.makedirs(path, exist_ok = True)

if config["remove_duplicated_reads"]:
    rm_dup_flag =  "--ignoreDup"
else:
    rm_dup_flag = ""

# Functions
def get_sra(wildcards):
            return [ my_files for my_files in glob.glob(f"data/00RawData/{wildcards.sample}.sra")]

# Set what rules to run locally
localrules: all #,
            #build_abundant_db

rule all:
    # IMPORTANT: output file fo all rule has to match the name specified in the output file
    # and not include suffixes that the command use might add to it.
    input:  o0 = expand("results/00reads/{s}.fastq.gz", s=samples),
            o1 = expand("results/01trim/{s}.trimmed.fastq.gz", s=samples),
            o2 = "results/05counts/read_counts.summary",
            o3 = expand("results/06fastqc_raw/{s}_fastqc.html", s=samples),        
            o4 = expand("results/06fastqc_trim/{s}.trimmed_fastqc.html", s=samples),
            o7 = "results/07multiqc/multiqc_done.flag",
            o9 = expand("results/05bigwig/{s}.bw", s=samples)

rule sra2fq:
    input: sra = get_sra
    output: "results/00reads/{sample}.fastq.gz"
    resources: 
        cpu_per_task = 1,
        partition = "quick",
        time = "3:00:00"
    shell:
        """
        fastq-dump --split-3 --gzip {input.sra} -O results/00reads/
        """

rule trimming:
    input:  fq1 = "results/00reads/{sample}.fastq.gz"
    output: fqU = "results/01trim/{sample}.trimmed.fastq.gz"
    params: 
            "ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=20 overwrite=t minlength=15"
    resources:
        cpus_per_task = 16,
        partition = "quick",
        time = "4:00:00"
    threads: 16
    log:    log1 = "results/01trim/{sample}.log",
            log2 = "results/01trim/{sample}.stats.log"
    benchmark:
            "benchmarks/trim/{sample}.tsv"
    shell:
        """
        bbduk.sh -Xmx1g threads={threads} \
            in={input.fq1} \
            out={output.fqU} \
            ref=data/00adapters/truseq.fa.gz \
            {params} stats={log.log2} 2> {log.log1}
        """

#rule build_abundant_db:
#    input: config["reference"]["abundant_rna_file"]
#    output: "data/00ref/abundant"
#    shell:
#        """
#        bowtie2-build {input} {output}
#        touch {output}
#        """

#rule rm_abundant_rnas:
#    input:  fq1P = "results/01trim/{sample}.1P.fastq.gz",
#            fq2P = "results/01trim/{sample}.2P.fastq.gz",
#            rnas_idx = "data/00ref/abundant"
#    output:
#        fq1F = "results/02abundant/{sample}.fastq.1.gz",
#        fq2F = "results/02abundant/{sample}.fastq.2.gz"
#    threads: 16
#    resources:
#        cpus_per_task = 16,
#        partition = "norm",
#        time = "24:00:00",
#        gres = "lscratch:20"
#    log: 
#        metrics = "results/02abundant/{sample}.metrics.txt",
#        logs = "results/02abundant/{sample}.log"
#    benchmark:
#            "benchmarks/rm_abundant_rnas/{sample}.tsv"
#    shell:
#        """
#        bowtie2 --threads {threads} -L 20 -x {input.rnas_idx} \
#         --met-file {log.metrics} \
#         --un-conc-gz results/02abundant/{wildcards.sample}.fastq.gz \
#         -1 {input.fq1P} -2 {input.fq2P} 2> {log.logs} 1> /dev/null
#        """ 
       
rule make_star_index:
    input: gen = f"{genome}", 
        ann = f"{annotation}" 
    output: db = "data/00ref/SA"
    resources: 
        mem_mb = 32000,
        cpus_per_task = 16,
        partition = "quick",
        time = "3:00:00",
        gres = "lscratch:20"
    threads: 16
    params: so = f"{starOverhang}",
            db_dir = "data/00ref"
    shell:
        """
            STAR --runMode genomeGenerate \
             --genomeDir {params.db_dir} \
             --genomeFastaFiles {input.gen} \
             --sjdbGTFfile {input.ann} \
             --sjdbOverhang {params.so}
            touch {output.db}  
        """

rule map_reads:
    input: fq1 = "results/01trim/{sample}.trimmed.fastq.gz",
           database = "data/00ref/SA"
    output: bam = "results/03map_reads/{sample}.Aligned.sortedByCoord.out.bam",
            log = "results/03map_reads/{sample}.Log.final.out"
    threads: 16
    resources:
        cpus_per_task = 16,
        partition = "norm",
        time = "14:00:00",
        mem_mb = 32000,
        gres = "lscratch:20"
    benchmark:
        "benchmarks/map_reads/{sample}.tsv"
    params: genome_dir = "data/00ref",
            prefix = "results/03map_reads/{sample}.",
    shell:
        """
        STAR --runMode alignReads \
                --runThreadN {threads} \
                --genomeDir {params.genome_dir} \
                --alignSJDBoverhangMin 1 \
                --alignSJoverhangMin 5 \
                --outFilterMismatchNmax 2 \
                --alignEndsType EndToEnd \
                --readFilesIn {input.fq1} \
                --readFilesCommand zcat --outFileNamePrefix {params.prefix} \
                --quantMode GeneCounts \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMattributes All
        touch {output.bam}
        touch {output.log}
        """

#rule remove_duplicates:
#    input: "results/03map_reads/{sample}.Aligned.sortedByCoord.out.bam"
#    output: "results/04dedup/{sample}.sorted.dedup.bam"
#    params: "READ_NAME_REGEX=null REMOVE_DUPLICATES=true"
#    log: "results/04dedup/{sample}.sorted.dedup.metrics.txt"
#    benchmark:
#        "benchmarks/remove_duplicates/{sample}.tsv"
#    resources:
#        cpus_per_task = 4,
#        mem_mb = 64000,
#        partition = "quick",
#        time = "4:00:00",
#        gres = "lscratch:20"
#    shell:
#        """
#        picard -Xmx32g MarkDuplicates \
#         I={input} \
#         O={output} \
#         M={log} \
#         {params}
#        samtools index {output}
#        """

rule remove_duplicates:
    input: "results/03map_reads/{sample}.Aligned.sortedByCoord.out.bam"
    output: "results/04dedup/{sample}.sorted.dedup.bam"
    params: "--READ_NAME_REGEX null --REMOVE_DUPLICATES false"
    log: "results/04dedup/{sample}.sorted.dedup.metrics.txt"
    threads: 24
    benchmark:
        "benchmarks/remove_duplicates/{sample}.tsv"
    resources:
        cpus_per_task = 4,
        mem_mb = 96000,
        partition = "norm",
        time = "4:00:00",
        gres = "lscratch:20"
    shell:
        """
        picard -Xmx16g -Xms16g -XX:ParallelGCThreads=5 MarkDuplicates \
         -I {input} \
         -O {output} \
         -M {log} \
         {params}
        samtools index {output}
        """


rule make_bigwig:
    input: "results/04dedup/{sample}.sorted.dedup.bam"
    output: "results/05bigwig/{sample}.bw"
    params: "--binSize 10 --normalizeUsing BPM" #  + "--filterRNAstrand [forward/reverse]" to plot strand-specific data
    threads: 8
    resources:
        cpus_per_task = 8,
        partition = "norm",
        time = "14:00:00"
    shell:
        """
        bamCoverage -p {threads} -b {input} -o {output} {params}
        """

rule run_correlation_analysis:
    input: expand("results/04dedup/{s}.sorted.dedup.bam", s=samples)
    output: "results/05correlation/multiBamSummary.results.npz"
    benchmark:
        "benchmarks/run_correlation_analysis/correlation.tsv"
    threads: 1
    resources:
        cpus_per_task = 1,
        partition = "norm",
        time = "14:00:00"
    shell: 
        """
        multiBamSummary bins --bamfiles {input} -o {output}
        plotCorrelation -in {output} --whatToPlot heatmap \
         --corMethod pearson -o results/05correlation/heatmap.png \
         --outFileCorMatrix outFileCorMatrix.txt
        """

rule counts:
    input: 
        genome = config["reference"]["genome_file"],
        gtf = config["reference"]["ensembl_gtf"],
        bam = expand("results/04dedup/{s}.sorted.dedup.bam", s=samples)
    output: counts = "results/05counts/read_counts",
            summary = "results/05counts/read_counts.summary"
    params: f"-t exon -g gene_id -O -s 2 -J -R BAM {rm_dup_flag} -M --fraction"  # Current params ignore multimappers and duplicated reads
                                                                   # -p = count fragments instead of individual reads
                                                                   # -M = include multi-mapping reads -O count reads mapping overlapping features
                                                                   # --fraction = multimapped reads will be caused as a fraction 
                                                                   #              instead of 1 (1/x where x = numb alignments reported for same read)
    benchmark:
        "benchmarks/counts/counts.tsv"
    threads: 16
    resources:
        cpus_per_task = 16,
        partition = "norm",
        time = "14:00:00",
        mem_mb = 32000
    shell:
        """
        featureCounts {params} -G {input.genome} -T {threads}\
         -a {input.gtf} \
         -o {output.counts} {input.bam}
        """

rule fastqc:
    input: raw = "results/00reads/{sample}.fastq.gz",
           trim = "results/01trim/{sample}.trimmed.fastq.gz"
    output: 
            o1 = "results/06fastqc_raw/{sample}_fastqc.html",
            o2 = "results/06fastqc_trim/{sample}.trimmed_fastqc.html",
            myzip1 = "results/06fastqc_raw/{sample}_fastqc.zip",
            myzip2 = "results/06fastqc_trim/{sample}.trimmed_fastqc.zip"
    threads: 8
    resources:
        cpus_per_task = 8,
        partition = "quick",
        time = "4:00:00",
        mem_mb = 4000
    params:
        "--quiet"
    shell:
        """
        fastqc {params} -t {threads} -o results/06fastqc_raw {input.raw}
        fastqc {params} -t {threads} -o results/06fastqc_trim {input.trim}
        touch {output.myzip1}
        touch {output.myzip2}
        """

rule multiqc:
    input: 
           i1 = expand("results/01trim/{sample}.log", sample = samples),
           i3 = expand("results/03map_reads/{sample}.Log.final.out", sample = samples),
           i4 = expand("results/04dedup/{sample}.sorted.dedup.metrics.txt", sample = samples),
           i5 = "results/05counts/read_counts.summary",
           i7 = expand("results/06fastqc_raw/{sample}_fastqc.zip", sample = samples),
           i8 = expand("results/06fastqc_trim/{sample}.trimmed_fastqc.zip", sample = samples)
    output: "results/07multiqc/multiqc_done.flag"
    resources:
        partition = "quick",
        time = "4:00:00",
        mem_mb = 4000                    
    shell:
        """
        multiqc -f -d -o results/07multiqc {input.i1} \
                {input.i7} \
                {input.i8} \
                {input.i3} \
                {input.i4} \
                {input.i5}
        touch {output}
        """


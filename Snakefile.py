import pandas as pd
from subprocess import run
from snakemake.utils import min_version
min_version("7.0")

# Uncomment for debugging
    # genomes = genomes[0]
    # kmers = kmers[:1]
# wildcard_constraints:
#     exp=experiments,
#     acc=acc_list

#Config
kmc_genome_count_jobs = 1 #Limits kmc_genome_counts rule to 1 job at a time

work_dir = "/zata/zippy/andrewsg/scratch"
genomes = ["GRCh38", "T2T"]
kmers = [str(_) for _ in range(50,101,5)]
metadata = pd.read_csv("Metadata.txt", sep="\t", header=0)
metadata = metadata[metadata["exp"] == "ENCSR590IHT"]

experiments = metadata.exp.unique().tolist()
acc_list = metadata.acc.tolist()
URLs = dict(zip(metadata.acc, metadata.url))
print(URLs)

# import rules
include: "./rules/genome.smk.py"
include: "./rules/fastq.smk.py"
include: "./rules/qc.smk.py"
include: "./rules/aln.smk.py"

rule all:
    input: 
        # expand(work_dir + "/genome/{genome}.fa", genome=genomes),
        # expand(work_dir + "/kmer/{genome}-k{kmer}.kmc_pre", genome=genomes, kmer=kmers),
        # expand(work_dir + "/kmer/{genome}-k{kmer}.bw", genome=genomes, kmer=kmers),
        # expand(work_dir + "/genome/{genome}.fa.1.bt2", genome=genomes),
        # expand(work_dir + "/fastq-raw/{acc}.fastq.gz", acc=acc_list)
        # expand(work_dir + "/fastq-merged/{exp}_{read}.fastq.gz", exp=experiments, read=[str(i) for i in range(1,3)]),
        # expand(work_dir + "/fastq-trimmed/{exp}_1.P.fastq.gz", exp=experiments),
        expand(work_dir + "/aln-raw/{genome}-{exp}.bam", genome=genomes, exp=experiments),
        # expand(work_dir + "/aligned/{genome}-{exp}.presorted.bam", genome=genomes, exp=experiments),
        # expand(work_dir + "/aligned/{genome}-{exp}.bam.bai", genome=genomes, exp=experiments),
        # expand(work_dir + "/aligned/{genome}-{exp}.filtered.bam", genome=genomes, exp=experiments),
        # expand(work_dir + "/aligned/{genome}-{exp}-k{kmer}.bam", genome=genomes, exp=experiments, kmer=kmers),
        # expand(work_dir + "/output/{genome}-{exp}-k{kmer}.postsort.bam", genome=genomes, exp=experiments, kmer=kmers)


    
#fastqc prescreen


rule samtools_index:
    input:
        presort = work_dir + "/aligned/{genome}-{exp}.presorted.bam",
    output:
        index = temp(work_dir + "/aligned/{genome}-{exp}.bam.bai"),
    threads: 32
    resources:
        mem_mb=240000,
        c=32,
        runtime=240,
        nodes=1,
        slurm_partition="4hours"
    log:
        work_dir + "/logs/samtools_index/{genome}-{exp}.log"
    singularity:
        "docker://clarity001/t2t-encode"
    shell:
        """
        (
        set -e
        set -o pipefail
        echo "Running samtools index"
        samtools index \
        -@{threads} \
        {input.presort} \
        {output.index}
        ) &> {log}
        """

rule samtools_filter:
    input:
        index = work_dir + "/aligned/{genome}-{exp}.bam.bai",
        presort = work_dir + "/aligned/{genome}-{exp}.presorted.bam",
    output:
        temp(work_dir + "/aligned/{genome}-{exp}.filtered.bam"),
    threads: 32
    resources:
        mem_mb=240000,
        c=32,
        runtime=240,
        nodes=1,
        slurm_partition="4hours"
    log:
        work_dir + "/logs/samtools-filter/{genome}-{exp}.log",
    singularity:
        "docker://clarity001/t2t-encode",
    shell:
        """
        (
        set -e
        set -o pipefail
        echo "Running samtools filter"
        samtools view -b -F 1804 -f 2 -q 2 \
        -@{threads} \
        {input.presort} > {output}
        ) &> {log}
        """


rule unique_kmer_filtering:
    input:
        bigWigs = expand(work_dir + "/kmer/{genome}-k{kmer}.bw", kmer=kmers)
        bam = work_dir + "/aln-filtered/{genome}-{exp}.filtered.bam",
    output:
        work_dir + "/aln-final/.bam"),
    threads: 1
    resources:
        mem_mb=4000,
        c=1,
        runtime=240,
        nodes=1,
        slurm_partition="4hours"
    log:
        work_dir + "/logs/unique_kmer_filtering/{genome}-{exp}-k{kmer}.log"
    singularity:
        "docker://clarity001/t2t-encode"
    params:
        template =  work_dir + "/BigWig/{genome}-k{kmer}-bw/{genome}-k*.bw",
        py_script = "/opt/T2T_Encode_Analysis/bin/filter_by_unique_kmers.py"
    shell:
        """
        (
        set -e
        set -o pipefail
        echo "Filtering unique kmers"
        python3 {params.py_script} {input.filtered} {params.template} {wildcards.kmer} {output}
        ) &> {log}
        """

rule samtools_postsort:
    input:
        work_dir + "/aligned/{genome}-{exp}-k{kmer}.bam",
    output:
        work_dir + "/output/{genome}-{exp}-k{kmer}.postsort.bam",
    threads:2
    resources:
        mem_mb=240000,
        c=32,
        runtime=240,
        nodes=1,
        slurm_partition="4hours"
    log:
        work_dir + "/logs/samtools_postsort/{genome}-{exp}-k{kmer}.log",
    singularity:
        "docker://clarity001/t2t-encode"
    shell:
        """
        (
        set -e
        set -o pipefail
        echo "Running samtools sort"
        samtools sort \
        -o {output} \
        {input} \
        --threads {threads}
        echo "Done"
        ) &> {log}
        """

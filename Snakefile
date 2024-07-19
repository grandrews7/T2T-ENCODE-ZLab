import pandas as pd
from subprocess import run
from snakemake.utils import min_version
min_version("6.0")

# Uncomment for debugging
    # genomes = genomes[0]
    # kmers = kmers[:1]
# wildcard_constraints:
#     exp=experiments,
#     acc=acc_list

#Config
kmc_genome_count_jobs = 1
work_dir = "/zata/zippy/ramirezc/scratch"
genomes = ["GRCh38", "T2T"]
kmers = [str(_) for _ in range(50,101,5)]
metadata = pd.read_csv("Metadata.txt", sep="\t", header=0, nrows=4)
experiments = metadata.exp.unique().tolist()
acc_list = metadata.acc.tolist()
URLs = dict(zip(metadata.acc, metadata.url))
print(URLs)

rule all:
    input: 
        expand(work_dir + "/kmer/{genome}-k{kmer}.bw", genome=genomes, kmer=kmers),
        expand(work_dir + "/genome/{genome}.fa.1.bt2", genome=genomes),
        expand(work_dir + "/fastq-raw/{acc}.fastq.gz", acc=acc_list),
        expand(work_dir + "/fastq-merged/{exp}_{read}.fastq.gz", exp=experiments, read=[str(i) for i in range(1,3)]),
        expand(work_dir + "/fastq-trimmed/{exp}_1.paired.fastq.gz", exp=experiments),
        expand(work_dir + "/aligned/{genome}-{exp}.bam", genome=genomes, exp=experiments),
        expand(work_dir + "/aligned/{genome}-{exp}.presorted.bam", genome=genomes, exp=experiments),
        expand(work_dir + "/aligned/{genome}-{exp}.bam.bai", genome=genomes, exp=experiments),
        expand(work_dir + "/aligned/{genome}-{exp}.filtered.bam", genome=genomes, exp=experiments),
        expand(work_dir + "/aligned/{genome}-{exp}-k{kmer}.bam", genome=genomes, exp=experiments, kmer=kmers),

ruleorder: get_fastq > merge_fastq

rule get_fastq:
    output: temp(work_dir + "/fastq-raw/{acc}.fastq.gz"),
    threads: 1
    log:
        work_dir + "/logs/get_fastq/{acc}.log",
    params: 
        url = lambda w: URLs[w.acc]
    shell:
        """
        (
        set -e
        set -o pipefail
        echo "Downloading fastq files"
        wget -q -O {output} {params.url}
        ) &> {log}
        """

def get_exp_fastq_list(wildcards):
    return([work_dir + "/fastq-raw/" + acc + ".fastq.gz" for acc in metadata[metadata["exp"] == wildcards.exp].acc.tolist()])

def get_R1(exp):
    return(metadata[(metadata["exp"] == exp) & (metadata["read"] == 1)].acc.tolist())

def get_R2(exp):
    return(metadata[(metadata["exp"] == exp) & (metadata["read"] == 2)].acc.tolist())
    
rule merge_fastq:
    input: get_exp_fastq_list
    output: 
        fq1 = work_dir + "/fastq-merged/{exp}_1.fastq.gz",
        fq2 = work_dir + "/fastq-merged/{exp}_2.fastq.gz",
    threads: 24
    # resources:
    #     mem_mb=240000,
    #     c=24,
    #     runtime=240,
    #     nodes=1,
    #     slurm_partition="4hours"
    log:
        work_dir + "/logs/merge_fastq/{exp}.log",
    singularity:
        "docker://andrewsg/t2t-encode",
    params:
        fastq_raw_dir = work_dir + "/fastq-raw",
        fastq_merged_dir = work_dir + "/fastq-merged",
        R1 = lambda w: get_R1(w.exp),
        R2 = lambda w: get_R2(w.exp),
    shell: 
        """
        (
        set -e
        set -o pipefail
        echo "Merging fastq files"
        for acc in {params.R1}; do
            zcat -1 {params.fastq_raw_dir}/$acc.fastq.gz >> {params.fastq_merged_dir}/{wildcards.exp}_1.fastq 
        done

        for acc in {params.R2}; do
            zcat -1 {params.fastq_raw_dir}/$acc.fastq.gz >> {params.fastq_merged_dir}/{wildcards.exp}_2.fastq
        done
        pigz -p {threads} {params.fastq_merged_dir}/{wildcards.exp}_1.fastq
        pigz -p {threads} {params.fastq_merged_dir}/{wildcards.exp}_2.fastq
        #touch {output.fq1}
        #touch {output.fq2}
        ) &> {log}
        """

#fastqc prescreen

rule trimmomatic:
    input: 
        fq1 = work_dir + "/fastq-merged/{exp}_1.fastq.gz",
        fq2 = work_dir + "/fastq-merged/{exp}_2.fastq.gz",
    output: 
        fq1_paired = temp(work_dir + "/fastq-trimmed/{exp}_1.paired.fastq.gz"),
        fq1_unpaired = temp(work_dir + "/fastq-trimmed/{exp}_1.unpaired.fastq.gz"),
        fq2_paired = temp(work_dir + "/fastq-trimmed/{exp}_2.paired.fastq.gz"),
        fq2_unpaired = temp(work_dir + "/fastq-trimmed/{exp}_2.unpaired.fastq.gz"),
    threads: 24
    # resources:
    #     mem_mb=240000,
    #     c=24,
    #     runtime=240,
    #     nodes=1,
    #     slurm_partition="4hours"
    log:
        work_dir + "/logs/trimmomatic/{exp}.log",
    singularity:
        "docker://andrewsg/t2t-encode",
    params:
        trimmomatic = "/opt/Trimmomatic-0.39/trimmomatic-0.39.jar",
        trim_method = "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36",
    shell:
        """
        (
        set -e
        set -o pipefail
        echo "Trimming fastq files"
        java -jar {params.trimmomatic} PE -phred33 -threads {threads} {input.fq1} {input.fq2} {output.fq1_paired} {output.fq1_unpaired} {output.fq2_paired} {output.fq2_unpaired} {params.trim_method}
        ) &> {log}
        """

rule get_genome:
    output:
        fa = work_dir + "/genome/{genome}.fa",
        sizes = work_dir + "/genome/{genome}.sizes.txt",
    threads: 1
    # resources:
    #     mem_mb=4000,
    #     c=1,
    #     runtime=240,
    #     nodes=1,
    #     slurm_partition="4hours"
    log:
        work_dir + "/logs/get_genome/{genome}.log"
    singularity:
        "docker://andrewsg/t2t-encode"
    shell:
        """
        (
        set -e
        set -o pipefail
        echo "Downloading genome fasta files"
        if [ {wildcards.genome} = 'GRCh38' ]
        then
            echo 'GRCh38'
            url='http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'
        else
            echo 'T2T'
            url='http://t2t.gi.ucsc.edu/chm13/hub/t2t-chm13-v1.0/genome/t2t-chm13-v1.0.fa.gz'
        fi
        wget -q -O {output.fa}.gz $url
        gzip -d -c {output.fa}.gz > {output.fa}
        rm {output.fa}.gz
        samtools faidx {output.fa}
        cut -f1,2 {output.fa}.fai > {output.sizes}
        ) &> {log}
        """

rule bowtie2_build:
    input:
        fa = work_dir + "/genome/{genome}.fa",
    output: multiext(work_dir + "/genome/{genome}.fa", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    threads: 24
    # resources:
    #     mem_mb=240000,
    #     c=24,
    #     runtime=240,
    #     nodes=1,
    #     slurm_partition="4hours"
    log:
        work_dir + "/logs/bowtie2_build/{genome}.log"
    singularity:
        "docker://andrewsg/t2t-encode"
    params: 
        prefix = work_dir
    shell:
        """
        (
        set -e
        set -o pipefail
        echo "Running bowtie2 build"
        cd {params.prefix}/genome
        bowtie2-build --threads {threads} {wildcards.genome}.fa {wildcards.genome}.fa
        ) &> {log}
        """

rule bowtie2_align:
    input:
        fq1_paired = work_dir + "/fastq-trimmed/{exp}_1.paired.fastq.gz",
        fq2_paired = work_dir + "/fastq-trimmed/{exp}_2.paired.fastq.gz",
    output:
        temp(work_dir + "/aligned/{genome}-{exp}.raw.bam"),
    threads: 24
    # resources:
    #     mem_mb=240000,
    #     c=24,
    #     runtime=240,
    #     nodes=1,
    #     slurm_partition="4hours"
    log:
        work_dir + "/logs/bowtie2_align/{genome}-{exp}.log"
    singularity:
        "docker://andrewsg/t2t-encode"
    shell:
        """
        (
        set -e
        set -o pipefail
        echo "Running bowtie2 align"
        bowtie2 --no-discordant --no-mixed --very-sensitive --no-unal --omit-sec-seq --xeq --reorder \
        --threads {threads} \
        -x {work_dir}/genome/{wildcards.genome} \
        -1 {input.fq1_paired} \
        -2 {input.fq2_paired} |\
        samtools view -@{threads} \
        -o {output}
        ) &> {log}
        """

rule samtools_presort:
    input:
        work_dir + "/aligned/{genome}-{exp}.raw.bam",
    output:
        temp(work_dir + "/aligned/{genome}-{exp}.presorted.bam")
    threads: 24
    resources:
        mem_mb=240000,
        c=24,
        runtime=240,
        nodes=1,
        slurm_partition="4hours",
    log:
        work_dir + "/logs/samtools_presort/{genome}-{exp}.log",
    singularity:
        "docker://andrewsg/t2t-encode",
    shell:
        """
        (
        set -e
        set -o pipefail
        echo "Running samtools sort"
        samtools sort -o {output} {input} --threads {threads}
        ) &> {log}
        """

rule samtools_index:
    input:
        presort = work_dir + "/aligned/{genome}-{exp}.presorted.bam",
    output:
        index = temp(work_dir + "/aligned/{genome}-{exp}.bam.bai"),
    threads: 24
    # resources:
    #     mem_mb=240000,
    #     c=24,
    #     runtime=240,
    #     nodes=1,
    #     slurm_partition="4hours"
    log:
        work_dir + "/logs/samtools_index/{genome}-{exp}.log"
    singularity:
        "docker://andrewsg/t2t-encode"
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
        work_dir + "/aligned/{genome}-{exp}.bam.bai",
    output:
        work_dir + "/aligned/{genome}-{exp}.filtered.bam"
    threads: 24
    # resources:
    #     mem_mb=240000,
    #     c=24,
    #     runtime=240,
    #     nodes=1,
    #     slurm_partition="4hours"
    log:
        work_dir + "/logs/samtools-filter/{genome}-{exp}.log",
    singularity:
        "docker://andrewsg/t2t-encode",
    params:
        presort = work_dir + "/aligned/{genome}-{exp}.presorted.bam",
    shell:
        """
        (
        set -e
        set -o pipefail
        echo "Running samtools filter"
        samtools view -b -F 1804 -f 2 -q 2 \
        -@{threads} \
        {params.presort} > {output}
        ) &> {log}
        """

rule kmc:
    input:
        fa = work_dir + "/genome/{genome}.fa",
    output:
        pre = work_dir + "/kmer/{genome}-k{kmer}.kmc_pre",
        shuf = work_dir + "/kmer/{genome}-k{kmer}.kmc_suf"
    threads: 24
    # resources:
    #     mem_mb=240000,
    #     c=24,
    #     runtime=240,
    #     nodes=1,
    #     slurm_partition="4hours"
    log:
        work_dir + "/logs/kmc/{genome}-k{kmer}.log"
    singularity:
        "docker://andrewsg/t2t-encode"
    params: prefix = work_dir
    shell:
        """
        (
        set -e
        set -o pipefail
        echo "Creating kmer suffix trees"
        rm -rf /tmp/{wildcards.genome}-k{wildcards.kmer}; mkdir /tmp/{wildcards.genome}-k{wildcards.kmer}
        kmc -k{wildcards.kmer} -m128 -fm -ci1 -t{threads} {input.fa} {params.prefix}/kmer/{wildcards.genome}-k{wildcards.kmer} /tmp/{wildcards.genome}-k{wildcards.kmer}
        ) &> {log}
        """

rule kmc_genome_counts: #ensure to run use these flags: --cores all --resources kmc_genome_count_jobs=1 
    input:
        fa = work_dir + "/genome/{genome}.fa",
        sizes = work_dir + "/genome/{genome}.sizes.txt",
        pre = work_dir + "/kmer/{genome}-k{kmer}.kmc_pre",
        shuf = work_dir + "/kmer/{genome}-k{kmer}.kmc_suf",
    output:
        wig = workDir + "kmer/{genome}-k{kmer}.wig",
    threads: 24
    resources:
    #     mem_mb=240000,
    #     c=16,
    #     runtime=240,
    #     nodes=1,
    #     slurm_partition="4hours"
        kmc_genome_count_jobs=1
    log:
        work_dir + "/logs/kmc_genome_counts/{genome}-k{kmer}.log"
    singularity:
        "docker://andrewsg/t2t-encode"
    params: 
        prefix = work_dir
    shell:
        """
        (
        set -e
        set -o pipefail
        echo 'Obtaining genome counts'
        kmc_genome_counts -t{threads} {params.prefix}/kmer/{wildcards.genome}-k{wildcards.kmer} {input.fa} > /tmp/{wildcards.genome}-k{wildcards.kmer}.wig
        """

rule wigToBigWig:
    input:
        wig = workDir + "kmer/{genome}-k{kmer}.wig",
        sizes = workDir + "genome/{genome}.sizes.txt",
    output:
        bw = workDir + "kmer/{genome}-k{kmer}.bw",
    threads: 24
    resources:
    #     mem_mb=240000,
    #     c=16,
    #     runtime=240,
    #     nodes=1,
    #     slurm_partition="4hours"
    log:
        work_dir + "/logs/wigToBigWig/{genome}-k{kmer}.log"
    singularity:
        "docker://andrewsg/t2t-encode"
    shell:
        """
        (
        set -e
        set -o pipefail
        echo 'Running wigToBigWig'
        wigToBigWig {input.wig} {input.sizes} {output.bw}
        ) &> {log}
        """

rule unique_kmer_filtering:
    input:
        bigWig = work_dir + "/kmer/{genome}-k{kmer}.bw",
        filtered = work_dir + "/aligned/{genome}-{exp}.filtered.bam",
    output:
        temp(work_dir + "/aligned/{genome}-{exp}-k{kmer}.bam"),
    threads: 1
    # resources:
    #     mem_mb=240000,
    #     c=24,
    #     runtime=240,
    #     nodes=1,
    #     slurm_partition="4hours"
    log:
        work_dir + "/logs/kmc_genome_counts/{genome}-{exp}-k{kmer}.log"
    singularity:
        "docker://andrewsg/t2t-encode"
    params:
        template = work_dir + "/kmer/{genome}-k*.bw",
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
    threads: 24
    # resources:
    #     mem_mb=240000,
    #     c=24,
    #     runtime=240,
    #     nodes=1,
    #     slurm_partition="4hours"
    log:
        work_dir + "/logs/samtools_postsort/{genome}-k{kmer}.log",
    singularity:
        "docker://andrewsg/t2t-encode"
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

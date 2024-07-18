import pandas as pd
from subprocess import run
workDir = "/zata/zippy/andrewsg/scratch"
genomes = ["GRCh38", "T2T"]
kmers = [str(_) for _ in range(50,101,5)]

# Uncomment for debugging
#genomes = genomes[:1]
#kmers = kmers[:1]

metadata = pd.read_csv("Metadata.txt", sep="\t", header=0, nrows=4)
experiments = metadata.exp.unique().tolist()
acc_list = metadata.acc.tolist()

# wildcard_constraints:
#     exp=experiments,
#     acc=acc_list

URLs = dict(zip(metadata.acc, metadata.url))
print(URLs)

rule all:
    input: 
        expand(workDir + "/kmer/{genome}-k{kmer}.bigWig", genome=genomes, kmer=kmers),
        expand(workDir + "/genome/{genome}.fa.1.bt2", genome=genomes),
        expand(workDir + "/fastq-raw/{acc}.fastq.gz", acc=acc_list),
        expand(workDir + "/fastq-merged/{exp}_{read}.fastq.gz", exp=experiments, read=[str(i) for i in range(1,3)])

# rule all:
#     input:
#         expand(workDir + "/genome/{genome}.fa.1.bt2", genome=genomes)


ruleorder: get_fastq > merge_fastq

rule get_fastq:
    output: temp(workDir + "/fastq-raw/{acc}.fastq.gz")
    threads: 1
    params: url = lambda w: URLs[w.acc]
    shell:
        """
        wget -q -O {output} {params.url}
        """

def get_exp_fastq_list(wildcards):
    return([workDir + "/fastq-raw/" + acc + ".fastq.gz" for acc in metadata[metadata["exp"] == wildcards.exp].acc.tolist()])

def get_R1(exp):
    return(metadata[(metadata["exp"] == exp) & (metadata["read"] == 1)].acc.tolist())

def get_R2(exp):
    return(metadata[(metadata["exp"] == exp) & (metadata["read"] == 2)].acc.tolist())
    
rule merge_fastq:
    input: get_exp_fastq_list
    output: 
        fq1 = workDir + "/fastq-merged/{exp}_1.fastq.gz",
        fq2 = workDir + "/fastq-merged/{exp}_2.fastq.gz"
    threads: 1
    params:
        fastq-rawDir = workDir + "/fastq-raw",
        fastq-mergedDir = workDir + "/fastq-merged",
        R1 = lambda w: get_R1(w.exp),
        R2 = lambda w: get_R2(w.exp)
    shell:
        """
        for acc in {params.R1}; do
            zcat {params.fastq-rawDir}/$acc.fastq.gz >> {params.fastq-mergedDir}/{wildcards.exp}_1.fastq
        done

        for acc in {params.R2}; do
            zcat {params.fastq-rawDir}/$acc.fastq.gz >> {params.fastq-mergedDir}/{wildcards.exp}_2.fastq
        done
        gzip {params.fastq-mergedDir}/{wildcards.exp}_1.fastq
        gzip {params.fastq-mergedDir}/{wildcards.exp}_2.fastq
        #touch {output.fq1}
        #touch {output.fq2}
        """

rule trimmomatic:
    input: 
        fq1 = workDir + "/fastq-merged/{exp}_1.fastq",
        fq2 = workDir + "/fastq-merged/{exp}_2.fastq"
    output: 
        fq1Paired = temp(workDir + "/fastq-trimmed/{exp}_1.paired.fastq")
        fq1Unpaired = temp(workDir + "/fastq-trimmed/{exp}_1.unpaired.fastq")
        fq1Paired = temp(workDir + "/fastq-trimmed/{exp}_2.paired.fastq")
        fq2Unpaired = temp(workDir + "/fastq-trimmed/{exp}_1.unpaired.fastq")
    threads: 24
    params:
        step = "string"
    shell:
        """
        java -jar trimmomatic-0.39 PE -phred64 -threads {threads} {input.fq1} {input.fq2} {output.fq1Paired} {output.fq1Unpaired} {output.fq2Paired} {output.fq2Unpaired} {params.step}
        """

rule get_genome:
    output:
        fa = workDir + "/genome/{genome}.fa",
        sizes = workDir + "/genome/{genome}.sizes.txt"

    singularity:
        "docker://andrewsg/t2t-encode"
    threads: 1
    resources:
        mem_mb=4000,
        c=1,
        runtime=240,
        nodes=1,
        slurm_partition="4hours"
    shell:
        """
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
        """

rule bowtie2_build:
    input:
        fa = workDir + "/genome/{genome}.fa"
    
    output: multiext(workDir + "/genome/{genome}.fa", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    params: prefix = workDir
    threads: 24
    resources:
        mem_mb=240000,
        c=24,
        runtime=240,
        nodes=1,
        slurm_partition="4hours"
    singularity:
        "docker://andrewsg/t2t-encode"
    shell:
        """
        cd {params.prefix}/genome
        bowtie2-build --threads {threads} {wildcards.genome}.fa {wildcards.genome}.fa
        """
    
rule kmc:
    input:
        fa = workDir + "/genome/{genome}.fa",
    output:
        pre = workDir + "/kmer/{genome}-k{kmer}.kmc_pre",
        shuf = workDir + "/kmer/{genome}-k{kmer}.kmc_suf"
    singularity:
        "docker://andrewsg/t2t-encode"
    threads: 24
    resources:
        mem_mb=240000,
        c=24,
        runtime=240,
        nodes=1,
        slurm_partition="4hours"
    params: prefix = workDir
    shell:
        """
        rm -rf /tmp/{wildcards.genome}-k{wildcards.kmer}; mkdir /tmp/{wildcards.genome}-k{wildcards.kmer}
        kmc -k{wildcards.kmer} -m128 -fm -ci1 -t{threads} {input.fa} {params.prefix}/kmer/{wildcards.genome}-k{wildcards.kmer} /tmp/{wildcards.genome}-k{wildcards.kmer}
        """

rule kmc_genome_counts:
    input:
        fa = workDir + "/genome/{genome}.fa",
        sizes = workDir + "/genome/{genome}.sizes.txt",
        pre = workDir + "/kmer/{genome}-k{kmer}.kmc_pre",
        shuf = workDir + "/kmer/{genome}-k{kmer}.kmc_suf"
    output:
        bigWig = workDir + "/kmer/{genome}-k{kmer}.bigWig",
    params: prefix = workDir
    threads: 16
    resources:
        mem_mb=240000,
        c=16,
        runtime=240,
        nodes=1,
        slurm_partition="4hours"
    singularity:
        "docker://andrewsg/t2t-encode"
    shell:
        """
        echo 'Obtaining genome counts'
        kmc_genome_counts -t{threads} {params.prefix}/kmer/{wildcards.genome}-k{wildcards.kmer} {input.fa} > /tmp/{wildcards.genome}-k{wildcards.kmer}.wig
        echo 'Running wigToBigWig'
        wigToBigWig /tmp/{wildcards.genome}-k{wildcards.kmer}.wig {input.sizes} {output}
        """


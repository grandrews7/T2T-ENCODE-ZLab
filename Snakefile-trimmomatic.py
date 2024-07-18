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
        expand(workDir + "/fastq1/{acc}.fastq.gz", acc=acc_list),
        expand(workDir + "/fastq2/{exp}_{read}.fastq.gz", exp=experiments, read=[str(i) for i in range(1,3)])

# rule all:
#     input:
#         expand(workDir + "/genome/{genome}.fa.1.bt2", genome=genomes)


ruleorder: get_fastq > merge_fastq

rule get_fastq:
    output: temp(workDir + "/fastq1/{acc}.fastq.gz")
    threads: 1
    params: url = lambda w: URLs[w.acc]
    shell:
        """
        wget -q -O {output} {params.url}
        """

def get_exp_fastq_list(wildcards):
    return([workDir + "/fastq1/" + acc + ".fastq.gz" for acc in metadata[metadata["exp"] == wildcards.exp].acc.tolist()])

def get_R1(exp):
    return(metadata[(metadata["exp"] == exp) & (metadata["read"] == 1)].acc.tolist())

def get_R2(exp):
    return(metadata[(metadata["exp"] == exp) & (metadata["read"] == 2)].acc.tolist())
    
rule merge_fastq:
    input: get_exp_fastq_list
    output: 
        fq1 = workDir + "/fastq2/{exp}_1.fastq.gz",
        fq2 = workDir + "/fastq2/{exp}_2.fastq.gz"
    threads: 1
    params:
        fastq1Dir = workDir + "/fastq1",
        fastq2Dir = workDir + "/fastq2",
        R1 = lambda w: get_R1(w.exp),
        R2 = lambda w: get_R2(w.exp)
    shell:
        """
        for acc in {params.R1}; do
            zcat {params.fastq1Dir}/$acc.fastq.gz >> {params.fastq2Dir}/{wildcards.exp}_1.fastq
        done

        for acc in {params.R2}; do
            zcat {params.fastq1Dir}/$acc.fastq.gz >> {params.fastq2Dir}/{wildcards.exp}_2.fastq
        done
        gzip {params.fastq2Dir}/{wildcards.exp}_1.fastq
        gzip {params.fastq2Dir}/{wildcards.exp}_2.fastq
        #touch {output.fq1}
        #touch {output.fq2}
        """
        
#add trimmomatic rule here

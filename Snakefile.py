import pandas as pd
from subprocess import run
from snakemake.utils import min_version
min_version("7.0")

work_dir = "/zata/zippy/andrewsg/scratch"
genomes = ["GRCh38", "T2T"]
kmers = [str(_) for _ in range(50,101,5)]

# Read in metadata
metadata = pd.read_csv("Metadata.txt", sep="\t", header=0)

# Experiments to test pipeline
# mixture of paired-ended and single-ended
# and replicated and unreplicated

test_experiments = ["ENCSR724YTA", "ENCSR841XDW", "ENCSR798NVH", "ENCSR617IFZ", "ENCSR822CEA"]
# test_experiments = ["ENCSR822CEA"]
metadata = metadata[metadata["exp"].isin(test_experiments)]
metadata_pe = metadata[metadata["run_type"] == "paired-ended"].reset_index(drop=True)
metadata_se = metadata[metadata["run_type"] == "single-ended"].reset_index(drop=True)


experiments = metadata.exp.unique().tolist()
print("There are {} experiments".format(len(experiments)))
pe_experiments = metadata[metadata["run_type"] == "paired-ended"].exp.unique().tolist() 
se_experiments = metadata[metadata["run_type"] == "single-ended"].exp.unique().tolist()
print("{} single-ended experiments".format(len(se_experiments)))
print("{} paired-ended experiments".format(len(pe_experiments)))

exp_bioreps = list(set([exp + "-" + str(biorep) for exp, biorep in zip(metadata.exp, metadata.biorep)]))
print("{} unique pairs of experiments and biological replicates".format(len(exp_bioreps)))
print(exp_bioreps)
exp_bioreps_pe = list(set([exp + "-" + str(biorep) for exp, biorep in zip(metadata_pe.exp, metadata_pe.biorep)]))
exp_bioreps_se = list(set([exp + "-" + str(biorep) for exp, biorep in zip(metadata_se.exp, metadata_se.biorep)]))
print("...{} are paired-ended".format(len(exp_bioreps_pe)))
print(exp_bioreps_pe)
print("...{} are single-ended".format(len(exp_bioreps_se)))
print(exp_bioreps_se)
acc_list = metadata.acc.tolist()
#acc_list = acc_list[:1]
URLs = dict(zip(metadata.acc, metadata.url))
MD5SUMs = dict(zip(metadata.acc, metadata.md5sum))

# import rules
# include: "./rules/genome.smk.py"
# include: "./rules/fastq.smk.py"
# include: "./rules/qc.smk.py"
# include: "./rules/aln.smk.py"

workdir: work_dir
rule all:
    input: 
        expand(work_dir + "/fastq-raw/{acc}.fastq.gz", acc = acc_list),
        [work_dir + "/fastq-merged-se/" + _ + ".fastq.gz" for _ in exp_bioreps_se],
        [work_dir + "/fastq-merged-pe/" + _ + "_1.fastq.gz" for _ in exp_bioreps_pe],
        [work_dir + "/fastq-merged-pe/" + _ + "_2.fastq.gz" for _ in exp_bioreps_pe],
        [work_dir + "/fastq-trimmed-pe/" + _ + "_1.P.fastq.gz" for _ in exp_bioreps_pe],
        [work_dir + "/fastq-trimmed-se/" + _ + ".fastq.gz" for _ in exp_bioreps_se]

# download individual fastq files
rule get_fastq:
    output: work_dir + "/fastq-raw/{acc}.fastq.gz",
    threads: 1
    log:
        work_dir + "/logs/get_fastq/{acc}.log"
    resources:
        mem_mb=4000,
        c=1,
        runtime=240,
        nodes=1,
        slurm_partition="4hours"
    singularity:
        "docker://andrewsg/t2t-encode",
    params: 
        url = lambda w: URLs[w.acc],
        md5sum = lambda w: MD5SUMs[w.acc]
    shell:
        """
        (
        echo Changing into working directory
        cd {workflow._workdir}
        echo "Downloading fastq files"
        checksum=""
        while [ "$checksum" != {params.md5sum} ]
        do
            wget -q -O {output} {params.url}
            checksum=$(md5sum {output} | awk '{{print $1; exit}}')
        done
        ) &> {log}
        """

def get_fastq_list_se(wildcards):
    return([work_dir + "/fastq-raw/" + acc + ".fastq.gz" for acc in metadata[(metadata["exp"] == wildcards.exp) & (metadata["biorep"] == int(wildcards.biorep))].acc.tolist()])
    
rule merge_fastq_se:
    input: get_fastq_list_se
    output: work_dir + "/fastq-merged-se/{exp}-{biorep}.fastq.gz" 
    threads: 32
    log:
        work_dir + "/logs/merge_fastq_se/{exp}-{biorep}.log"
    shell:
        """
        echo Changing into working directory
        cd {workflow._workdir}
        echo {input}
        pigz -cd -p {threads} {input} | pigz -c -p {threads} > {output} 
        """

def get_fastq_list_pe(wildcards):
    tmp_df = metadata[(metadata["exp"] == wildcards.exp) & (metadata["biorep"] == int(wildcards.biorep)) & (metadata["read"] == int(wildcards.read))]
    return([work_dir + "/fastq-raw/" + acc + ".fastq.gz" for acc in tmp_df.sort_values("read_count").acc.tolist()])
    
rule merge_fastq_pe:
    input: get_fastq_list_pe
    output: work_dir + "/fastq-merged-pe/{exp}-{biorep}_{read}.fastq.gz"
    threads: 32
    log:
        work_dir + "/logs/merge_fastq_pe/{exp}-{biorep}_{read}.log"
    shell:
        """
        echo Changing into working directory
        cd {workflow._workdir}
        echo {input}
        pigz -cd -p {threads} {input} | pigz -c -p {threads} > {output} 
        """

rule trim_pe:
    input: 
        R1 = work_dir + "/fastq-merged-pe/{exp}-{biorep}_1.fastq.gz",
        R2 = work_dir + "/fastq-merged-pe/{exp}-{biorep}_2.fastq.gz"
    output:
        fq1_paired = work_dir + "/fastq-trimmed-pe/{exp}-{biorep}_1.P.fastq.gz",
        fq1_unpaired = work_dir + "/fastq-trimmed-pe/{exp}-{biorep}_1.U.fastq.gz",
        fq2_paired = work_dir + "/fastq-trimmed-pe/{exp}-{biorep}_2.P.fastq.gz",
        fq2_unpaired = work_dir + "/fastq-trimmed-pe/{exp}-{biorep}_2.U.fastq.gz",
    threads: 64
    resources:
        mem_mb=90000,
        c=24,
        runtime=240,
        nodes=1,
        slurm_partition="4hours"
    log:
        work_dir + "/logs/trim_pe/{exp}-{biorep}.log",
    singularity:
        "docker://andrewsg/t2t-encode",
    params:
        trimmomatic = "/opt/Trimmomatic-0.39/trimmomatic-0.39.jar",
        trim_method = "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36",
    shell:
        """
        (
        echo Changing into working directory
        cd {workflow._workdir}
        echo "Trimming fastq files"
        java -jar {params.trimmomatic} PE -phred33 -threads {threads} {input.R1} {input.R2} {output.fq1_paired} {output.fq1_unpaired} {output.fq2_paired} {output.fq2_unpaired} {params.trim_method}
        ) &> {log}
        """

rule trim_se:
    input: work_dir + "/fastq-merged-se/{exp}-{biorep}.fastq.gz",
    output: work_dir + "/fastq-trimmed-se/{exp}-{biorep}.fastq.gz",
    threads: 64
    resources:
        mem_mb=90000,
        c=24,
        runtime=240,
        nodes=1,
        slurm_partition="4hours"
    log:
        work_dir + "/logs/trim_se/{exp}-{biorep}.log",
    singularity:
        "docker://andrewsg/t2t-encode",
    params:
        trimmomatic = "/opt/Trimmomatic-0.39/trimmomatic-0.39.jar",
        trim_method = "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36",
    shell:
        """
        (
        echo Changing into working directory
        cd {workflow._workdir}
        echo "Trimming fastq files"
        java -jar {params.trimmomatic} SE -phred33 -threads {threads} {input} {output} {params.trim_method}
        ) &> {log}
        """



























        
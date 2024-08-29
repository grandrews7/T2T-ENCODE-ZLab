import pandas as pd
from subprocess import run
from snakemake.utils import min_version
min_version("8.0")
import os
import re
import sys

snakeDir = os.getcwd()
scriptDir = snakeDir + "/scripts"
print(scriptDir)
workdir: "/zata/zippy/andrewsg/scratch"
workDir = os.getcwd()
basedir = "/data/zusers/andrewsg/TF-STR/T2T-ENCODE-Pipeline"
genomes = ["GRCh38", "T2T"]
#genomes = ["GRCh38"]
kmers = [str(_) for _ in range(30,101,5)]
# Read in metadata
metadata = pd.read_csv(basedir + "/Metadata.txt", sep="\t", header=0)
ZNF_TFs = metadata[metadata["target"].str.contains("ZNF")].target.unique().tolist()
print("There are {} ZNF TFs".format(len(ZNF_TFs)))

trt_experiments = metadata[metadata["target"].isin(ZNF_TFs)].exp.unique().tolist()
ctrl_experiments = metadata[metadata["exp"].isin(trt_experiments)].control.unique().tolist()

experiments = trt_experiments + ctrl_experiments


# qc = pd.read_csv(basedir + "/QC.txt", sep="\t", header=0)
# experiments = list(set([x.split("-")[1] for x in qc[qc["aln-raw"] == 0]["Genome-Experiment-Biorep"].tolist()]))
# experiments = ["ENCSR063FUV", "ENCSR934WOF"]
metadata = metadata[metadata["exp"].isin(experiments)]

metadata_pe = metadata[metadata["run_type"] == "paired-ended"].reset_index(drop=True)
metadata_se = metadata[metadata["run_type"] == "single-ended"].reset_index(drop=True)
print("Paired-ended metadata")
print(metadata_pe.head())
print("Single-ended metadata")
print(metadata_se.head())

CONTROLS = dict(zip(metadata[metadata["target"] != "Control"].exp, metadata[metadata["target"] != "Control"].control))

print("There are {} experiments".format(len(experiments)))
pe_experiments = metadata[metadata["run_type"] == "paired-ended"].exp.unique().tolist() 
se_experiments = metadata[metadata["run_type"] == "single-ended"].exp.unique().tolist()
print("{} single-ended experiments".format(len(se_experiments)))
print("{} paired-ended experiments".format(len(pe_experiments)))

exp_bioreps = list(set([exp + "-" + str(biorep) for exp, biorep in zip(metadata.exp, metadata.biorep)]))

print("{} unique pairs of experiments and biological replicates".format(len(exp_bioreps)))
exp_bioreps_pe = list(set([exp + "-" + str(biorep) for exp, biorep in zip(metadata_pe.exp, metadata_pe.biorep)]))
exp_bioreps_se = list(set([exp + "-" + str(biorep) for exp, biorep in zip(metadata_se.exp, metadata_se.biorep)]))
print("...{} are paired-ended".format(len(exp_bioreps_pe)))
print("...{} are single-ended".format(len(exp_bioreps_se)))
acc_list = metadata.acc.tolist()
#acc_list = acc_list[:1]
URLs = dict(zip(metadata.acc, metadata.url))
MD5SUMs = dict(zip(metadata.acc, metadata.md5sum))
RUN_TYPES = dict(zip(metadata.exp, metadata.run_type))
# import rules
include: "./rules/kmer.smk.py"
include: "./rules/genome.smk.py"
include: "./rules/fastq.smk.py"
include: "./rules/aln.smk.py"
include: "./rules/qc.smk.py"
rule all:
    input:
        # expand(workDir + "/genome/{genome}.fa.1.bt2", genome=genomes),
        # expand(workDir + "/fastq-raw/{acc}.fastq.gz", acc = acc_list),
        # expand(workDir + "/kmer/{genome}-k{kmer}.bw", genome=genomes, kmer=kmers),
        # [workDir + "/fastq-merged-se/" + _ + ".fastq.gz" for _ in exp_bioreps_se],
        # [workDir + "/fastq-merged-pe/" + _ + "_1.fastq.gz" for _ in exp_bioreps_pe],
        # [workDir + "/fastq-merged-pe/" + _ + "_2.fastq.gz" for _ in exp_bioreps_pe],
        # [workDir + "/fastq-trimmed-se/" + _ + ".fastq.gz" for _ in exp_bioreps_se],
        # [workDir + "/fastq-trimmed-pe/" + _ + "_1.P.fastq.gz" for _ in exp_bioreps_pe],
        # [workDir + "/aln-raw-se/" + g + "-" + _ + ".bam" for _ in exp_bioreps_se for g in genomes],
        #[workDir + "/aln-filtered/" + g + "-" + _ + ".bam" for _ in exp_bioreps for g in genomes],
        # [workDir + "/aln-final/" + g + "-" + _ + ".bam" for _ in exp_bioreps for g in genomes],
        [workDir + "/qc/aln-final/" + g + "-" + _ + ".flagstat" for _ in exp_bioreps for g in genomes],
        # expand(workDir + "/peaks/{genome}-{exp}.narrowPeak", exp = trt_experiments)



    
    


    




























        
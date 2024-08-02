rule get_chain_file:
    output:
        chain1 = workDir + "/genome/grch38-chm13v2.chain",
        chain2 = workDir + "/genome/chm13v2-grch38.chain",
    params:
        url1 = "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/grch38-chm13v2.chain",
        url2 = "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-grch38.chain"
    log:
        workDir + "/logs/get-chain-file.log"
    threads: 1
    resources:
        mem_mb = 16000,
        runtime = 240,
        c = 1,
        nodes = 1,
        slurm_partition = "4hours"
    singularity:
        "docker://andrewsg/t2t-encode"
    shell:
        """
        (
        echo "Downloading chain files"
        wget -q -O {output.chain1} {params.url1} 
        wget -q -O {output.chain2} {params.url2}
        ) &> {log}
        """

def get_chain_file(wildcards):
    if wildcards.genome == "GRCh38":
        return workDir + "/genome/grch38-chm13v2.chain"
    else:
        return workDir + "/genome/chm13v2-grch38.chain",

rule liftOver:
    input:
        sumpeaks = workDir + "/macs2/{genome}-{exp}_peaks.narrowPeak",
        chainfile = get_chain_file
    output:
        bed = workDir + "/liftOver/{genome}-{exp}.{lift_type}.bed",
    params:
        min_ratio = 0.95
    log:
        workDir + "/logs/liftOver/{genome}-{exp}.log"
    threads: 1
    resources:
        mem_mb = 16000,
        runtime = 240,
        c = 1,
        nodes = 1,
        slurm_partition = "4hours"
    singularity:
        "docker://andrewsg/t2t-encode"
    shell: # bedPlus=6 treat only the first 6 columns in your input file as a BED file
        """
        (
        echo "Lifting over"
        liftOver \
        -bedPlus=6 \
        {input.sumpeaks} \
        {input.chainfile} \
        {workDir}/liftOver/{genome}-{exp}.shared.bed \
        {workDir}/liftOver/{genome}-{exp}.only.bed 
        ) &> {log}
        """
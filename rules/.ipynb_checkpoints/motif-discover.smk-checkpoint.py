rule get_fasta:
    input:
        in_fa = workDir + "/genome/{genome}.fa",
        bed = workDir + "/liftOver/{genome}-{exp}.{lift_type}.bed",
    output:
        out_fa = workDir + "/STREME/{genome}-{exp}.{lift_type}.fa",
    log:
         workDir + "/logs/STREME/get_fasta/{genome}-{exp}.{lift_type}.log"
    threads: 1
    resources:
        mem_mb = 16000,
        runtime = 240,
        c = 1,
        nodes = 1,
        slurm_partition = "4hours"
    singularity:
        "docker://clarity001/t2t-encode"
    shell: 
        """
        (
        echo "Generating fasta file"
        bedtools getfasta -fi {input.fa} -bed {input.bed}
        ) &> {log}
        """

rule_motif_discover
    input:
        workDir + "/STREME/{genome}-{exp}.{lift_type}.fa",
    output:
        workDir + "/STREME/{genome}-{exp}.{lift_type}.txt"
    log:
        workDir + "/logs/STREME/motif_discover/{genome}-{exp}.{lift_type}.log"
    threads: 1
    resources:
        mem_mb = 16000,
        runtime = 240,
        c = 1,
        nodes = 1,
        slurm_partition = "4hours"
    singularity:
        "docker://clarity001/t2t-encode"
    shell: 
        """
        streme --nmotifs 8 --order 2 --text -p {input} > {output} 2> {log}
        """
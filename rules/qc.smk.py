 rule flagstat:
    input:
        workDir + "/aln-final/{genome}-{exp}-{biorep}.bam"
    
    output:
        workDir + "/qc/aln-final/{genome}-{exp}-{biorep}.flagstat",

    threads: 8
    resources:
        mem_mb=32000,
        c=8,
        runtime=720,
        nodes=1,
        slurm_partition="12hours"
    log:
        workDir + "/logs/flagstat/{genome}-{exp}-{biorep}.log"
    singularity:
        "docker://andrewsg/t2t-encode"
    shell:
        """
        (
        samtools flagstat {input} > {output}
        ) &> {log}
        """
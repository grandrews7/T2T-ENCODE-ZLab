def get_genome_size(wildcards):
    if wildcards.genome == "GRCh38":
        return "2913022398"
    else:
        return "3117292070"

rule macs2_callpeak:
    input:
        treatment = workDir + "/aln-final/{genome}-{exp}.bam",
        control = lambda w: workDir + "/aln-final/{genome}-" + CONTROLS[w.exp] + ".bam"
    output:
        peaks = workDir + "/macs2/{genome}-{exp}_peaks.narrowPeak",
    params:
        name = "{genome}-{exp}",
        outdir = workDir + "/macs2",
        genome_size = get_genome_size,
        run_type = lambda w: RUN_TYPES[w.exp]
    log:
        workDir + "/logs/macs2/{genome}-{exp}.log"
    threads: 32
    resources:
        mem_mb = 32000,
        runtime = 240,
        c = 32,
        nodes = 1,
        slurm_partition = "4hours"
    singularity:
        "docker://clarity001/t2t-encode"
    shell:
        """
        (
        echo "Calling peaks"
        if [ {params.run_type} == "paired-ended" ]; then
            echo "Sorting paired-ended BAM by name"
            samtools sort -@ {threads} -o {wildcards.genome}-{wildcards.exp}.namesort.bam {input.treatment}
            macs3 callpeak -t {wildcards.genome}-{wildcards.exp}.namesort.bam -c {input.control} \
                -f BAMPE -g {params.genome_size} -n {params.name} -B -q 0.01 \
                --outdir {params.outdir}
            rm {wildcards.genome}-{wildcards.exp}.namesort.bam
        else
            macs3 callpeak -t {input.treatment} -c {input.control} \
                -f BAM -g {params.genome_size} -n {params.name} -B -q 0.01 \
                --outdir {params.outdir}
        fi
        ) &> {log}
        """
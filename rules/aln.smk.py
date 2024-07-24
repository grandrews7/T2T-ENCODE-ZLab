rule bowtie2_align_se:
    input:
        work_dir + "/fastq-trimmed-se/{exp}-{biorep}.fastq.gz"
    
    output:
        work_dir + "/aln-raw-se/{genome}-{exp}-{biorep}.bam"

    threads: 24
    resources:
        mem_mb=90000,
        c=24,
        runtime=240,
        nodes=1,
        slurm_partition="12hours"
    log:
        work_dir + "/logs/bowtie2_align_se/{genome}-{exp}-{biorep}.log"
    singularity:
        "docker://clarity001/t2t-encode"
    shell:
        """
        (
        echo "Running bowtie2 align"
        bowtie2 --no-discordant --no-mixed --very-sensitive --no-unal --omit-sec-seq --xeq --reorder \
        --threads {threads} \
        -x {work_dir}/genome/{wildcards.genome}.fa \
        -1 {input[0]} \
        -2 {input[1]} |\
        samtools view -@{threads} \
        -o {output}
        ) &> {log}
        """

    
rule bowtie2_align_pe:
    input:
        work_dir + "/fastq-trimmed/{exp}-{biorep}_1.P.fastq.gz", 
        work_dir + "/fastq-trimmed/{exp}-{biorep}_2.P.fastq.gz"

    output:
        work_dir + "/aln-raw/{genome}-{exp}-{biorep}.bam",
    threads: 64
    resources:
        mem_mb=240000,
        c=64,
        runtime=240,
        nodes=1,
        slurm_partition="4hours"
    log:
        work_dir + "/logs/bowtie2_align/{genome}-{exp}-{biorep}.log"
    singularity:
        "docker://clarity001/t2t-encode"
    shell:
        """
        (
        echo "Running bowtie2 align"
        bowtie2 --no-discordant --no-mixed --very-sensitive --no-unal --omit-sec-seq --xeq --reorder \
        --threads {threads} \
        -x {work_dir}/genome/{wildcards.genome}.fa \
        -1 {input[0]} \
        -2 {input[1]} |\
        samtools view -@{threads} \
        -o {output}
        ) &> {log}
        """

rule samtools_sort:
    input:
        work_dir + "/aln-raw/{genome}-{exp}-{biorep}.bam",
    output:
        work_dir + "/aln-coordsorted/{genome}-{exp}-{biorep}.bam",
    threads: 24
    resources:
        mem_mb=90000,
        c=24,
        runtime=240,
        nodes=1,
        slurm_partition="4hours",
    log:
        work_dir + "/logs/sort/{genome}-{exp}-{biorep}.log",
    singularity:
        "docker://clarity001/t2t-encode",
    shell:
        """
        (
        echo "Running samtools sort"
        samtools sort -@ {threads} -o {output} {input}
        ) &> {log}
        """
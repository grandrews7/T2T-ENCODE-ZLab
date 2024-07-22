rule bowtie2_align:
    input:
        fq1_paired = work_dir + "/fastq-trimmed/{exp}_1.P.fastq.gz",
        fq2_paired = work_dir + "/fastq-trimmed/{exp}_2.P.fastq.gz",
    output:
        work_dir + "/aln-raw/{genome}-{exp}.bam",
    threads: 64
    resources:
        mem_mb=240000,
        c=64,
        runtime=240,
        nodes=1,
        slurm_partition="4hours"
    log:
        work_dir + "/logs/bowtie2_align/{genome}-{exp}.log"
    singularity:
        "docker://clarity001/t2t-encode"
    shell:
        """
        (
        echo "Running bowtie2 align"
        bowtie2 --no-discordant --no-mixed --very-sensitive --no-unal --omit-sec-seq --xeq --reorder \
        --threads {threads} \
        -x {work_dir}/genome/{wildcards.genome}.fa \
        -1 {input.fq1_paired} \
        -2 {input.fq2_paired} |\
        samtools view -@{threads} \
        -o {output}
        ) &> {log}
        """

rule samtools_sort:
    input:
        work_dir + "/aln-raw/{genome}-{exp}.bam",
    output:
        temp(work_dir + "/aln-raw/{genome}-{exp}.presorted.bam"),
    threads: 24
    resources:
        mem_mb=90000,
        c=24,
        runtime=240,
        nodes=1,
        slurm_partition="4hours",
    log:
        work_dir + "/logs/samtools_presort/{genome}-{exp}.log",
    singularity:
        "docker://clarity001/t2t-encode",
    shell:
        """
        (
        echo "Running samtools sort"
        samtools sort -@ {threads} -o {output} {input}
        ) &> {log}
        """
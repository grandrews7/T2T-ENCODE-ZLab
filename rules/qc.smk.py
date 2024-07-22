rule trim:
    input: 
        fq1 = work_dir + "/fastq-merged/{exp}_1.fastq",
        fq2 = work_dir + "/fastq-merged/{exp}_2.fastq",
    output: 
        fq1_paired = work_dir + "/fastq-trimmed/{exp}_1.P.fastq.gz",
        fq1_unpaired = work_dir + "/fastq-trimmed/{exp}_1.U.fastq.gz",
        fq2_paired = work_dir + "/fastq-trimmed/{exp}_2.P.fastq.gz",
        fq2_unpaired = work_dir + "/fastq-trimmed/{exp}_2.U.fastq.gz",
    threads: 24
    resources:
        mem_mb=90000,
        c=24,
        runtime=240,
        nodes=1,
        slurm_partition="4hours"
    log:
        work_dir + "/logs/trimmomatic/{exp}.log",
    singularity:
        "docker://clarity001/t2t-encode",
    params:
        trimmomatic = "/opt/Trimmomatic-0.39/trimmomatic-0.39.jar",
        trim_method = "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36",
    shell:
        """
        (
        echo "Trimming fastq files"
        java -jar {params.trimmomatic} PE -phred33 -threads {threads} {input.fq1} {input.fq2} {output.fq1_paired} {output.fq1_unpaired} {output.fq2_paired} {output.fq2_unpaired} {params.trim_method}
        ) &> {log}
        """

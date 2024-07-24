def get_fastq_list(wildcards):
    print(wildcards.exp, wildcards.biorep)
    tmp = metadata[(metadata["exp"] == wildcards.exp) & (metadata["biorep"] == int(wildcards.biorep))].reset_index(drop=True)
    tmp = tmp.sort_values("read", ascending=True)
    print(tmp)
    return([work_dir + "/fastq-raw/" + acc + ".fastq.gz" for acc in tmp.acc.tolist()])


rule trim_pe:
    input: 
        get_fastq_list
    output: 
        fq1_paired = work_dir + "/fastq-trimmed/{exp}-{biorep}_1.P.fastq.gz",
        fq1_unpaired = work_dir + "/fastq-trimmed/{exp}-{biorep}_1.U.fastq.gz",
        fq2_paired = work_dir + "/fastq-trimmed/{exp}-{biorep}_2.P.fastq.gz",
        fq2_unpaired = work_dir + "/fastq-trimmed/{exp}-{biorep}_2.U.fastq.gz",
    threads: 24
    resources:
        mem_mb=90000,
        c=24,
        runtime=240,
        nodes=1,
        slurm_partition="4hours"
    log:
        work_dir + "/logs/trimmomatic/{exp}-{biorep}.log",
    singularity:
        "docker://clarity001/t2t-encode",
    params:
        trimmomatic = "/opt/Trimmomatic-0.39/trimmomatic-0.39.jar",
        trim_method = "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36",
    shell:
        """
        (
        echo "Trimming fastq files"
        java -jar {params.trimmomatic} PE -phred33 -threads {threads} {input[0]} {input[1]} {output.fq1_paired} {output.fq1_unpaired} {output.fq2_paired} {output.fq2_unpaired} {params.trim_method}
        ) &> {log}
        """

# download individual fastq files
rule get_fastq:
    output: workDir + "/fastq-raw/{acc}.fastq.gz",
    threads: 1
    log:
        workDir + "/logs/get_fastq/{acc}.log"
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
    return([workDir + "/fastq-raw/" + acc + ".fastq.gz" for acc in metadata[(metadata["exp"] == wildcards.exp) & (metadata["biorep"] == int(wildcards.biorep))].acc.tolist()])
    
rule merge_fastq_se:
    input: get_fastq_list_se
    output: workDir + "/fastq-merged-se/{exp}-{biorep}.fastq.gz" 
    threads: 8
    resources:
        mem_mb=32000,
        c=8,
        runtime=240,
        nodes=1,
        slurm_partition="4hours"
    log:
        workDir + "/logs/merge_fastq_se/{exp}-{biorep}.log"
    params:
        inp_len = lambda wildcards, input: len(input)
    shell:
        """
        (
        echo {input}
        echo {params.inp_len}
        if [[ {params.inp_len} -eq 1 ]]; then
            mv {input} {output}
        else
            pigz -cd -p {threads} {input} | pigz -c -p {threads} > {output} 
        fi
        ) &> {log}
        """

def get_fastq_list_pe(wildcards):
    tmp_df = metadata[(metadata["exp"] == wildcards.exp) & (metadata["biorep"] == int(wildcards.biorep)) & (metadata["read"] == int(wildcards.read))]
    return([workDir + "/fastq-raw/" + acc + ".fastq.gz" for acc in tmp_df.sort_values("read_count").acc.tolist()])
    
rule merge_fastq_pe:
    input: get_fastq_list_pe
    output: workDir + "/fastq-merged-pe/{exp}-{biorep}_{read}.fastq.gz"
    threads: 8
    resources:
        mem_mb=32000,
        c=8,
        runtime=240,
        nodes=1,
        slurm_partition="4hours"
    log:
        workDir + "/logs/merge_fastq_pe/{exp}-{biorep}_{read}.log"
    params:
        inp_len =lambda wildcards, input: len(input)
    shell:
        """
        (
        echo {input}
        echo {params.inp_len}
        if [[ {params.inp_len} -eq 1 ]]; then
            mv {input} {output}
        else
            pigz -cd -p {threads} {input} | pigz -c -p {threads} > {output} 
        fi
        ) &> {log}
        """

rule trim_pe:
    input: 
        R1 = workDir + "/fastq-merged-pe/{exp}-{biorep}_1.fastq.gz",
        R2 = workDir + "/fastq-merged-pe/{exp}-{biorep}_2.fastq.gz"
    output:
        fq1_paired = workDir + "/fastq-trimmed-pe/{exp}-{biorep}_1.P.fastq.gz",
        fq1_unpaired = workDir + "/fastq-trimmed-pe/{exp}-{biorep}_1.U.fastq.gz",
        fq2_paired = workDir + "/fastq-trimmed-pe/{exp}-{biorep}_2.P.fastq.gz",
        fq2_unpaired = workDir + "/fastq-trimmed-pe/{exp}-{biorep}_2.U.fastq.gz"
    threads: 8
    resources:
        mem_mb=32000,
        c=8,
        runtime=240,
        nodes=1,
        slurm_partition="4hours"
    log:
        workDir + "/logs/trim_pe/{exp}-{biorep}.log",
    singularity:
        "docker://andrewsg/t2t-encode",
    params:
        trimmomatic = "/opt/Trimmomatic-0.39/trimmomatic-0.39.jar",
        trim_method = "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25",
    shell:
        """
        (
        echo "Trimming fastq files"
        java -jar {params.trimmomatic} PE -phred33 -threads {threads} {input.R1} {input.R2} {output.fq1_paired} {output.fq1_unpaired} {output.fq2_paired} {output.fq2_unpaired} {params.trim_method}
        ) &> {log}
        """

rule trim_se:
    input: workDir + "/fastq-merged-se/{exp}-{biorep}.fastq.gz",
    output: workDir + "/fastq-trimmed-se/{exp}-{biorep}.fastq.gz",
    threads: 8
    resources:
        mem_mb=32000,
        c=24,
        runtime=240,
        nodes=1,
        slurm_partition="4hours" 
    log:
        workDir + "/logs/trim_se/{exp}-{biorep}.log",
    singularity:
        "docker://andrewsg/t2t-encode",
    params:
        trimmomatic = "/opt/Trimmomatic-0.39/trimmomatic-0.39.jar",
        trim_method = "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25",
    shell:
        """
        (
        echo "Trimming fastq files"
        java -jar {params.trimmomatic} SE -phred33 -threads {threads} {input} {output} {params.trim_method}
        ) &> {log}
        """
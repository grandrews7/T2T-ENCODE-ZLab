def get_exp_fastq_list(wildcards):
    return([work_dir + "/fastq-raw/" + acc + ".fastq.gz" for acc in metadata[metadata["exp"] == wildcards.exp].acc.tolist()])

def get_reads(exp, read):
    return(metadata[(metadata["exp"] == exp) & (metadata["read"] == int(read))].acc.tolist())

    
rule get_fastq:
    output: temp(work_dir + "/fastq-raw/{acc}.fastq.gz"),
    threads: 1
    resources:
        mem_mb=4000,
        c=1,
        runtime=240,
        nodes=1,
        slurm_partition="4hours"
    log:
        work_dir + "/logs/get_fastq/{acc}.log",
    singularity:
        "docker://clarity001/t2t-encode",
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

rule merge_fastq:
    input: get_exp_fastq_list
    output: 
        temp(work_dir + "/fastq-merged/{exp}_{read}.fastq")
    threads: 8
    resources:
        mem_mb=16000,
        c=4,
        runtime=240,
        nodes=1,
        slurm_partition="4hours"
    log:
        work_dir + "/logs/merge_fastq/{exp}-{read}.log",
    singularity:
        "docker://clarity001/t2t-encode",
    params:
        fastq_raw_dir = work_dir + "/fastq-raw",
        fastq_merged_dir = work_dir + "/fastq-merged",
        reads = lambda w: get_reads(w.exp, w.read)
    shell: 
        """
        (
        echo "Merging fastq files"
        echo {wildcards}
        echo {params}
        for acc in {params.reads}; do
            echo $acc
            pigz -cd -p {threads} {params.fastq_raw_dir}/$acc.fastq.gz >> {params.fastq_merged_dir}/{wildcards.exp}_{wildcards.read}.fastq
        done
        ) &> {log}
        """
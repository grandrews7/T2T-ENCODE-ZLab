rule get_fastq:
    output: work_dir + "/fastq-raw/{acc}.fastq.gz",
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
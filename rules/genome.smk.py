rule get_genome:
    output:
        fa = work_dir + "/genome/{genome}.fa",
        sizes = work_dir + "/genome/{genome}.sizes.txt",
    threads: 1
    resources:
        mem_mb=4000,
        c=1,
        runtime=240,
        nodes=1,
        slurm_partition="4hours"
    log:
        work_dir + "/logs/get_genome/{genome}.log"
    singularity:
        "docker://clarity001/t2t-encode"
    shell:
        """
        (
        echo "Downloading genome fasta files"
        if [ {wildcards.genome} = 'GRCh38' ]
        then
            echo 'GRCh38'
            url='http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'
        else
            echo 'T2T'
            url='http://t2t.gi.ucsc.edu/chm13/hub/t2t-chm13-v1.0/genome/t2t-chm13-v1.0.fa.gz'
        fi
        wget -q -O {output.fa}.gz $url
        gzip -d -c {output.fa}.gz > {output.fa}
        rm {output.fa}.gz
        samtools faidx {output.fa}
        cut -f1,2 {output.fa}.fai > {output.sizes}
        ) &> {log}
        """

rule bowtie2_build:
    input:
        fa = work_dir + "/genome/{genome}.fa",
    output:
        multiext(work_dir + "/genome/{genome}.fa", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    threads: 32
    resources:
        mem_mb=240000,
        c=32,
        runtime=240,
        nodes=1,
        slurm_partition="4hours"
    log:
        work_dir + "/logs/bowtie2_build/{genome}.log"
    singularity:
        "docker://clarity001/t2t-encode"
    params: 
        prefix = work_dir
    shell:
        """
        (
        set -e
        set -o pipefail
        echo "Running bowtie2 build"
        cd {params.prefix}/genome
        bowtie2-build --threads {threads} {wildcards.genome}.fa {wildcards.genome}.fa
        ) &> {log}
        """
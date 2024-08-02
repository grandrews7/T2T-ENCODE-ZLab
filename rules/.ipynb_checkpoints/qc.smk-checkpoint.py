rule readcounts:
    output:
        qc = workDir + "/qc/raw_fastq/raw-fastq-readcounts.txt"
    threads: 1  
    resources:
        mem_mb=4000, 
        slurm_partition="30mins"
    run:
        import pandas as pd
        
        readcounts = metadata.groupby("exp")["read_count"].sum().to_dict()
        with open(output.qc, 'w') as f:
            for experiment, readcount_sum in readcounts.items():
                f.write(f"Raw FASTQ readcounts for {experiment}: {readcount_sum}\n")

rule qc_trim:
    input:
        trimmed_fastq = lambda w: expand(workDir + "/fastq-trimmed-se/{exp}-{biorep}.fastq.gz", exp=w.exp, biorep=metadata[(metadata.exp == w.exp)].biorep.unique())
                if w.exp in se_experiments else
                expand(workDir + "/fastq-trimmed-pe/{exp}-{biorep}_1.P.fastq.gz", exp=w.exp, biorep=metadata[(metadata.exp == w.exp)].biorep.unique()) +
                expand(workDir + "/fastq-trimmed-pe/{exp}-{biorep}_2.P.fastq.gz", exp=w.exp, biorep=metadata[(metadata.exp == w.exp)].biorep.unique())
    output:
        qc = workDir + "/qc/trim/{exp}.txt"
    threads: 32
    resources:
        mem_mb=32000,
        runtime=240,
        slurm_partition="4hours"
    singularity:
        "docker://andrewsg/t2t-encode"
    shell:
        """
        total_reads=0
        for fq in {input.trimmed_fastq}; do
            count=$(pigz -cd -p {threads} $fq | echo $((`wc -l`/4)))
            total_reads=$((total_reads + count))
        done
        echo "Trimmed FASTQ read count: $total_reads" > {output.qc}
        """

rule qc_bt2_align:
    input:
        bam = lambda w: expand(workDir + "/aln-raw-se/{genome}-{exp}-{biorep}.bam", genome=w.genome, exp=w.exp, biorep=metadata[(metadata.exp == w.exp)].biorep.unique())
                if w.exp in se_experiments else
                expand(workDir + "/aln-raw-pe/{genome}-{exp}-{biorep}.bam", genome=w.genome, exp=w.exp, biorep=metadata[(metadata.exp == w.exp)].biorep.unique())
    output:
        qc = workDir + "/qc/bt2_align/{genome}-{exp}.txt"
    threads: 4
    resources:
        mem_mb=8000,
        runtime=60,
        slurm_partition="30mins"
    singularity:
        "docker://andrewsg/t2t-encode"
    shell:
        """
        total_reads=0
        for bam in {input.bam}; do
            count=$(samtools flagstat $bam | grep 'in total' | awk '{{print $1}}')
            total_reads=$((total_reads + count))
        done
        echo "BT2 aligned read count: $total_reads" > {output.qc}
        """

rule qc_samtools_filter:
    input:
        bam = workDir + "/aln-filtered/{genome}-{exp}.bam"
    output:
        qc = workDir + "/qc/samtools_filter/{genome}-{exp}.txt"
    threads: 4
    resources:
        mem_mb=8000,
        runtime=60,
        slurm_partition="30mins"
    singularity:
        "docker://andrewsg/t2t-encode"
    shell:
        """
        samtools flagstat {input.bam} | grep 'in total' | awk '{{print "Samtools filtered read count:", $1}}' > {output.qc}
        """

rule qc_kmer_filter:
    input:
        bam = workDir + "/aln-final/{genome}-{exp}.bam"
    output:
        qc = workDir + "/qc/kmer_filter/{genome}-{exp}.txt"
    threads: 4
    resources:
        mem_mb=8000,
        runtime=60,
        slurm_partition="30mins"
    singularity:
        "docker://andrewsg/t2t-encode"
    shell:
        """
        samtools flagstat {input.bam} | grep 'in total' | awk '{{print "Kmer filtered read count:", $1}}' > {output.qc}
        """

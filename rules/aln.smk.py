rule bowtie2_align_se:
    input:
        workDir + "/fastq-trimmed-se/{exp}-{biorep}.fastq.gz",
    
    output:
        workDir + "/aln-raw-se/{genome}-{exp}-{biorep}.bam",

    threads: 24
    resources:
        mem_mb=32000,
        c=8,
        runtime=720,
        nodes=1,
        slurm_partition="12hours"
    log:
        workDir + "/logs/bowtie2_align_se/{genome}-{exp}-{biorep}.log"
    singularity:
        "docker://andrewsg/t2t-encode"
    shell:
        """
        (
        echo "Running bowtie2 align"
        /opt/bowtie2/bowtie2 --no-discordant --no-mixed --very-sensitive --no-unal --omit-sec-seq --xeq --reorder \
        --threads {threads} \
        -x /zata/zippy/andrewsg/scratch/genome/{wildcards.genome}.fa \
        --rg-id {wildcards.exp}-{wildcards.biorep} \
        --rg SM:{wildcards.exp}-{wildcards.biorep} \
        -U {input} | \
        samtools view -@{threads} \
        -o {output}
        ) &> {log}
        """

    
rule bowtie2_align_pe:
    input:
        workDir + "/fastq-trimmed-pe/{exp}-{biorep}_1.P.fastq.gz", 
        workDir + "/fastq-trimmed-pe/{exp}-{biorep}_2.P.fastq.gz",

    output:
        workDir + "/aln-raw-pe/{genome}-{exp}-{biorep}.bam",
    threads: 8
    resources:
        mem_mb=32000,
        c=8,
        runtime=720,
        nodes=1,
        slurm_partition="12hours"   
    log:
        workDir + "/logs/bowtie2_align_pe/{genome}-{exp}-{biorep}.log"
    singularity:
        "docker://andrewsg/t2t-encode"
    shell:
        """
        (
        echo "Running bowtie2 align"
        /opt/bowtie2/bowtie2 --no-discordant --no-mixed --very-sensitive --no-unal --omit-sec-seq --xeq --reorder \
        --threads {threads} \
        -x /zata/zippy/andrewsg/scratch/genome/{wildcards.genome}.fa \
        --rg-id {wildcards.exp}-{wildcards.biorep} \
        --rg SM:{wildcards.exp}-{wildcards.biorep} \
        -1 {input[0]} \
        -2 {input[1]} | \
        samtools view -@{threads} \
        -o {output}
        ) &> {log}
        """

def get_bams(wildcards):
    tmp_df = metadata[(metadata["exp"] == wildcards.exp)].reset_index(drop=True)
    read_type = tmp_df.run_type.tolist()[0]
    if read_type == "paired-ended":
        dir = workDir + "/aln-raw-pe"
    else:
        dir = workDir + "/aln-raw-se"

    bams = []
    bioreps = tmp_df.biorep.unique().tolist()
    for b in bioreps:
        bams.append(dir + "/" + wildcards.genome + "-" + wildcards.exp + "-" + str(b) + ".bam")
    print(bams)
    return(bams)

def get_flag(w):
    run_type = RUN_TYPES[w.exp]
    if run_type == "paired-ended":
        return(["1804", "2"])
    else:
        return(["1796", "0"])
        
    
    
rule filter_bam:
    input: get_bams
    output: workDir + "/aln-filtered/{genome}-{exp}.bam"
    threads: 32
    resources:
        mem_mb=32000,
        c=8,
        runtime=720,
        nodes=1,
        slurm_partition="12hours"  
    log:
        workDir + "/logs/filter_bam/{genome}-{exp}.log"
    singularity:
        "docker://andrewsg/t2t-encode"
    params:
            flag = lambda w: get_flag
    shell:
        """
        (
        echo "Merging BAM files"
        samtools merge -f -@ {threads} -o {wildcards.genome}-{wildcards.exp}.merge.bam {input}
        echo "Sorting by coordinate"
        samtools sort -@ {threads} -o {wildcards.genome}-{wildcards.exp}.sort.bam {wildcards.genome}-{wildcards.exp}.merge.bam
        echo "Indexing BAM file"
        samtools index -@ {threads} {wildcards.genome}-{wildcards.exp}.sort.bam
        

        java -jar /opt/picard/build/libs/picard.jar MarkDuplicates \
            INPUT={wildcards.genome}-{wildcards.exp}.sort.bam \
            OUTPUT={wildcards.genome}-{wildcards.exp}.markdup.bam \
            METRICS_FILE={wildcards.genome}-{wildcards.exp}.markdup.txt \
            ASSUME_SORT_ORDER=coordinate

        echo "Filtering BAM file (-f {params.flag[0]}; -F {params.flag[1]})"
        samtools view -b -F {params.flag[0]} -f {params.flag[1]} -q 2 \
        -@{threads} \
        -o  {output} \
        {wildcards.genome}-{wildcards.exp}.markdup.bam
        # rm {wildcards.genome}-{wildcards.exp}.merge.bam
        # rm {wildcards.genome}-{wildcards.exp}.sort.bam
        # rm {wildcards.genome}-{wildcards.exp}.markdup.bam
        ) &> {log}
        """


# Create dictionary of experiments and their read type (single-ended or paired-end)

rule kmer_filter:
    input: 
        bam = workDir + "/aln-filtered/{genome}-{exp}.bam",
        kmer_bigWigs = [workDir + "/kmer/{genome}-k" + _ + ".bw" for _ in kmers]
    output:
        workDir + "/aln-final/{genome}-{exp}.bam"
    threads: 32
    resources:
        mem_mb=32000,
        c=8,
        runtime=720,
        nodes=1,
        slurm_partition="12hours"   
    log:
        workDir + "/logs/kmer_filter/{genome}-{exp}.log"
    singularity:
        "docker://andrewsg/t2t-encode"
    params:
        run_type = lambda w: RUN_TYPES[w.exp],
        script = scriptDir + "/filter_by_unique_kmers.py",
        kmer_bigWigs = lambda w, input: ",".join(input.kmer_bigWigs)
    shell:
        """
        if [ {params.run_type} == "paired-ended" ]; then
            echo "Sorting paired-ended BAM by name"
            samtools sort -n -@ {threads} -o {wildcards.genome}-{wildcards.exp}.namesort.bam {input.bam}
            python3 {params.script} -b {wildcards.genome}-{wildcards.exp}.namesort.bam \
            -k {params.kmer_bigWigs} \
            -r {params.run_type} \
            -o {output}
        else
            echo "single-ended"
            python3 {params.script} -b {input.bam} \
            -k {params.kmer_bigWigs} \
            -r {params.run_type} \
            -o {output}
        fi
        """
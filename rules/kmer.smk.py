rule kmc:
    input:
        fa = work_dir + "/genome/{genome}.fa",
    output:
        pre = work_dir + "/kmer/{genome}-k{kmer}.kmc_pre",
        shuf = work_dir + "/kmer/{genome}-k{kmer}.kmc_suf"
    threads: 4
    resources:
        mem_mb=240000,
        c=4,
        runtime=240,
        nodes=1,
        slurm_partition="4hours"
    log:
        work_dir + "/logs/kmc/{genome}-k{kmer}.log"
    singularity:
        "docker://clarity001/t2t-encode"
    params: prefix = work_dir
    shell:
        """
        (
        echo "Creating kmer suffix trees"
        rm -rf /tmp/{wildcards.genome}-k{wildcards.kmer}; mkdir /tmp/{wildcards.genome}-k{wildcards.kmer}
        kmc -k{wildcards.kmer} -m128 -fm -ci1 -t{threads} {input.fa} {params.prefix}/kmer/{wildcards.genome}-k{wildcards.kmer} /tmp/{wildcards.genome}-k{wildcards.kmer}
        ) &> {log}
        """

rule kmc_genome_counts: #ensure to use flags --cores all --resources kmc_genome_count_jobs=1 when running snakemake
    input:
        fa = work_dir + "/genome/{genome}.fa",
        sizes = work_dir + "/genome/{genome}.sizes.txt",
        pre = work_dir + "/kmer/{genome}-k{kmer}.kmc_pre",
        shuf = work_dir + "/kmer/{genome}-k{kmer}.kmc_suf",
    output:
        wig = temp(work_dir + "/kmer/{genome}-k{kmer}.wig"),
    threads: 24
    resources:
        mem_mb=300000,
        c=24,
        runtime=240,
        nodes=1,
        slurm_partition="4hours",
        kmc_genome_count_jobs=1,
    log:
        work_dir + "/logs/kmc_genome_counts/{genome}-k{kmer}.log",
    singularity:
        "docker://clarity001/t2t-encode"
    params: 
        prefix = work_dir
    shell:
        """
        (
        echo 'Obtaining genome counts'
        kmc_genome_counts -t{threads} {params.prefix}/kmer/{wildcards.genome}-k{wildcards.kmer} {input.fa} > {output.wig}
        ) &> {log}
        """

rule wigToBigWig:
    input:
        wig = work_dir + "/kmer/{genome}-k{kmer}.wig",
        sizes = work_dir + "/genome/{genome}.sizes.txt",
    output:
        bw = work_dir + "/BigWig/{genome}-k{kmer}-bw/{genome}-k{kmer}.bw",
    threads: 1
    resources:
        mem_mb=4000,
        c=1,
        runtime=240,
        nodes=1,
        slurm_partition="4hours",
    log:
        work_dir + "/logs/wigToBigWig/{genome}-k{kmer}.log"
    singularity:
        "docker://clarity001/t2t-encode"
    shell:
        """
        (
        echo 'Running wigToBigWig'
        wigToBigWig {input.wig} {input.sizes} {output.bw}
        ) &> {log}
        """
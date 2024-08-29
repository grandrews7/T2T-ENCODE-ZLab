rule kmc:
    input:
        fa = workDir + "/genome/{genome}.fa",
    output:
        pre = workDir + "/kmer/{genome}-k{kmer}.kmc_pre",
        shuf = workDir + "/kmer/{genome}-k{kmer}.kmc_suf"
    threads: 8
    resources:
        mem_mb=240000,
        c=8,
        runtime=720,
        nodes=1,
        slurm_partition="5days"
    log:
        workDir + "/logs/kmc/{genome}-k{kmer}.log"
    singularity:
        "docker://andrewsg/t2t-encode"
    params: prefix = workDir
    shell:
        """
        (
        echo "Creating kmer suffix trees"
        rm -rf /tmp/{wildcards.genome}-k{wildcards.kmer}; mkdir /tmp/{wildcards.genome}-k{wildcards.kmer}
        kmc -k{wildcards.kmer} -m64 -fm -ci1 -t{threads} {input.fa} {params.prefix}/kmer/{wildcards.genome}-k{wildcards.kmer} /tmp/{wildcards.genome}-k{wildcards.kmer}
        ) &> {log}
        """

rule kmc_genome_counts:
    input:
        fa = workDir + "/genome/{genome}.fa",
        sizes = workDir + "/genome/{genome}.sizes.txt",
        pre = workDir + "/kmer/{genome}-k{kmer}.kmc_pre",
        shuf = workDir + "/kmer/{genome}-k{kmer}.kmc_suf",
    output:
        wig = workDir + "/kmer/{genome}-k{kmer}.wig",
    threads: 8
    resources:
        mem_mb=240000,
        c=8,
        runtime=720,
        nodes=1,
        slurm_partition="5days",
    log:
        workDir + "/logs/kmc_genome_counts/{genome}-k{kmer}.log",
    singularity:
        "docker://andrewsg/t2t-encode"
    params: 
        prefix = workDir
    shell:
        """
        (
        echo 'Obtaining genome counts'
        kmc_genome_counts -t{threads} {params.prefix}/kmer/{wildcards.genome}-k{wildcards.kmer} {input.fa} > {output.wig}
        ) &> {log}
        """

rule wigToBigWig:
    input:
        wig = workDir + "/kmer/{genome}-k{kmer}.wig",
        sizes = workDir + "/genome/{genome}.sizes.txt",
    output:
        bw = workDir + "/kmer/{genome}-k{kmer}.bw",
    threads: 8
    resources:
        mem_mb=240000,
        c=8,
        runtime=720,
        nodes=1,
        slurm_partition="5days",
    log:
        workDir + "/logs/wigToBigWig/{genome}-k{kmer}.log"
    singularity:
        "docker://andrewsg/t2t-encode"
    shell:
        """
        (
        echo 'Running wigToBigWig'
        wigToBigWig {input.wig} {input.sizes} {output.bw}
        ) &> {log}
        """
    
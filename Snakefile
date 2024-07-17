workDir = "/zata/zippy/andrewsg/scratch/"
genomes = ["GRCh38", "T2T"]
kmers = [str(_) for _ in range(50,101,5)]

# Uncomment for debugging
genomes = genomes[:1]
kmers = kmers[:1]

print(genomes, kmers)

# rule all:
#     input: 
#         expand(workDir + "kmer/{genome}-k{kmer}.bigWig", genome=genomes, kmer=kmers)

rule all:
    input:
        expand(workDir + "genome/{genome}.fa", genome=genomes)

rule download:
    output:
        fa = workDir + "genome/{genome}.fa",
        sizes = workDir + "genome/{genome}.sizes.txt"

    singularity:
        "docker://andrewsg/t2t-encode"
    threads: 1
    resources:
        mem_mb=4000,
        c=1,
        runtime=240,
        nodes=1,
        slurm_partition="4hours"
    shell:
        """
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
        """


rule kmc:
    input:
        fa = workDir + "genome/{genome}.fa",
    output:
        pre = workDir + "kmer/{genome}-k{kmer}.kmc_pre",
        shuf = workDir + "kmer/{genome}-k{kmer}.kmc_suf"
    singularity:
        "docker://andrewsg/t2t-encode"
    threads: 24
    resources:
        mem_mb=90000,
        c=24,
        runtime=240,
        nodes=1,
        slurm_partition="4hours"
    params: prefix = workDir
    shell:
        """
        rm -rf /tmp/{wildcards.genome}-k{wildcards.kmer}; mkdir /tmp/{wildcards.genome}-k{wildcards.kmer}
        kmc -k{wildcards.kmer} -m64 -fm -ci1 -t{threads} {input.fa} {params.prefix}kmer/{wildcards.genome}-k{wildcards.kmer} /tmp/{wildcards.genome}-k{wildcards.kmer}
        """

rule kmc_genome_counts:
    input:
        fa = workDir + "genome/{genome}.fa",
        sizes = workDir + "genome/{genome}.sizes.txt",
        pre = workDir + "kmer/{genome}-k{kmer}.kmc_pre",
        shuf = workDir + "kmer/{genome}-k{kmer}.kmc_suf"
    output:
        bigWig = workDir + "kmer/{genome}-k{kmer}.bigWig",
    params: prefix = workDir
    threads: 16
    resources:
        mem_mb=240000,
        c=16,
        runtime=240,
        nodes=1,
        slurm_partition="4hours"
    singularity:
        "docker://andrewsg/t2t-encode"
    shell:
        """
        echo 'Obtaining genome counts'
        kmc_genome_counts -t{threads} {params.prefix}kmer/{wildcards.genome}-k{wildcards.kmer} {input.fa} > /tmp/{wildcards.genome}-k{wildcards.kmer}.wig
        echo 'Running wigToBigWig'
        wigToBigWig /tmp/{wildcards.genome}-k{wildcards.kmer}.wig {input.sizes} {output}
        """


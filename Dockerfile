FROM ubuntu:22.04
MAINTAINER Greg Andrews <gregory.andrews@umassmed.edu>
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ="America/New_York"

RUN apt-get update && \
    apt-get install -y \
    build-essential \
    python3 python3-dev python3-pip \
    wget \
    curl \
    git \
    bedtools \
    bowtie2 \
    # for samtools 
    zlib1g-dev \
    libncurses-dev \
    libbz2-dev \
    liblzma-dev

RUN for i in wigToBigWig liftOver bigBedToBed bedToBigBed; do \
     wget -q http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/$i -O /bin/$i ; chmod +x /bin/$i ; done

RUN cd /opt && \
    git clone https://github.com/msauria/KMC.git && \
    cd KMC && \
    git checkout kmer_mapping && \
    make && \ 
    chmod a+rx bin/kmc && \
    chmod a+rx bin/kmc_genome_counts
    

RUN cd /usr/bin && \
     wget https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2 && \
     tar -vxjf samtools-1.20.tar.bz2 && \
     cd samtools-1.20 && \
     ./configure --prefix=/usr/bin/samtools && \
     make && \
     make install

RUN apt-get update && \
    apt-get install -y \
    openjdk-17-jdk openjdk-17-jre
    
RUN cd /opt && \
    git clone https://github.com/broadinstitute/picard.git && \
    cd picard/ && \
    ./gradlew shadowJar

# Add Trimmomatic, FastQC
# RUN git clone https://github.com/usadellab/Trimmomatic && \
# cd Trimmomatic
     
ENV PATH="/usr/bin/samtools/bin:/opt/KMC:/opt/KMC/bin:/opt/picard/build/libs:${PATH}"




    

FROM ubuntu:22.04
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ="America/New_York"

RUN apt-get update && \
    apt-get install -y \
    build-essential \
    python3 python3-dev python3-pip \
    perl \
    wget \
    curl \
    git \
    zip \
    pigz \
    bedtools \
    # for samtools 
    zlib1g-dev \
    libncurses-dev \
    libbz2-dev \
    liblzma-dev \
    default-jre \
    openjdk-17-jdk openjdk-17-jre

RUN pip install --upgrade pip && \
    pip install pysam pyBigWig numpy macs3 deeptools

RUN for i in wigToBigWig liftOver bigBedToBed bedToBigBed; do \
     wget -q http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/$i -O /bin/$i ; chmod +x /bin/$i ; done

RUN cd /opt && \
    git clone https://github.com/msauria/KMC.git && \
    cd KMC && \
    git checkout kmer_mapping && \
    make && \ 
    chmod a+rx bin/kmc && \
    chmod a+rx bin/kmc_genome_counts

ENV SAMTOOLS_VER="1.20"
RUN cd /usr/bin && \
     wget https://github.com/samtools/samtools/releases/download/"$SAMTOOLS_VER"/samtools-"$SAMTOOLS_VER".tar.bz2 && \
     tar -vxjf samtools-"$SAMTOOLS_VER".tar.bz2 && \
     cd samtools-"$SAMTOOLS_VER" && \
     ./configure --prefix=/usr/bin/samtools && \
     make && \
     make install
    
RUN cd /opt && \
    git clone https://github.com/broadinstitute/picard.git && \
    cd picard/ && \
    ./gradlew shadowJar

RUN cd /opt && \
    wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && \
	unzip Trimmomatic-0.39.zip && \
	echo "alias trimmomatic='java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar'" >> ~/.bash_aliases && \
    echo "source ~/.bash_aliases" >> ~/.bashrc && \
    rm Trimmomatic-0.39.zip && \
    chmod +x /opt/Trimmomatic-0.39/trimmomatic-0.39.jar

RUN cd /opt && \
    wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && \
    unzip fastqc_v0.12.1.zip && \
    cd FastQC && \
    chmod +x fastqc

ENV BOWTIE2_VER="2.5.4"

RUN wget -q -O bowtie2.zip http://sourceforge.net/projects/bowtie-bio/files/bowtie2/"$BOWTIE2_VER"/bowtie2-"$BOWTIE2_VER"-linux-x86_64.zip/download; \ 
	unzip bowtie2.zip -d /opt/; \ 
	ln -s /opt/bowtie2-"$BOWTIE2_VER"-linux-x86_64/ /opt/bowtie2; \ 
	rm bowtie2.zip 

RUN mkdir /opt/meme
ADD http://meme-suite.org/meme-software/5.4.1/meme-5.4.1.tar.gz /opt/meme
WORKDIR /opt/meme/
RUN tar zxvf meme-5.4.1.tar.gz && rm -fv meme-5.4.1.tar.gz
RUN cd /opt/meme/meme-5.4.1 && \
     ./configure --prefix=/opt  --enable-build-libxml2 --enable-build-libxslt  && \
     make && \
     make install && \
     rm -rfv /opt/meme

ENV PATH="/opt/libexec/meme-5.4.1:/opt/bin:${PATH}"
    
ENV PATH="/usr/bin/samtools/bin:/opt/KMC:/opt/KMC/bin:/opt/picard/build/libs:/opt/FastQC:/opt/Trimmomatic-0.39:${PATH}"
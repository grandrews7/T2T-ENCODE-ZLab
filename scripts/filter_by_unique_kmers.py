#!/usr/bin/env python3
print("Importing necessary libraries")
import sys
import pysam
import pyBigWig
import numpy
import argparse
import re
import os

def list_str(values):
    return values.split(',')
    
def arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bam", 
                        help="BAM file to filter",
                        required=True)
    
    parser.add_argument("-k", "--kmer_bigWigs",
                        help="Comma separated list of kmer bigWigs",
                        required = True,
                        type=list_str)
    
    parser.add_argument("-r", "--run_type", 
                        help = "Run type",
                        choices = ["single-ended", "paired-ended"],
                        required=True)
    
    parser.add_argument("-o", "--output", 
                        help = "Output BAM", 
                        required=True)

    
    return(parser)
    
def main():
    args = arg_parser().parse_args()
    bamFile = args.bam
    kmer_bigWigs = args.kmer_bigWigs
    output = args.output
    run_type = args.run_type
    sizes = [int(re.search("[k](\d){2,3}", _).group(0).lstrip("k")) for _ in kmer_bigWigs]
    sizes = numpy.array(sizes, dtype=numpy.int32)
    sizes2 = sizes[1:-1]
    DBs = {}
    for i, size in enumerate(sizes):
        DBs[size] = pyBigWig.open(kmer_bigWigs[i])

    bam = pysam.AlignmentFile(bamFile, 'rb')
    output = pysam.AlignmentFile(output, 'wb', template=bam)
    
    if run_type == "paired-ended":
        print("Filtering paired-ended BAM by unique kmers")
        bam_iter = bam.__iter__()
        try:
            read1 = next(bam_iter)
            while True:
                read2 = next(bam_iter)
                chrom = read1.reference_name
                if (valid_read(read1, chrom, DBs, sizes, sizes2)
                    or valid_read(read2, chrom, DBs, sizes, sizes2)):
                    output.write(read1)
                    output.write(read2)
                read1 = next(bam_iter)
        except StopIteration:
            pass
        finally:
            del bam_iter
        output.close()
        bam.close()
    else:
        print("Filtering single-ended BAM by unique kmers")
        for read in bam:
            chrom = read.reference_name
            if (valid_read(read, chrom, DBs, sizes, sizes2)):
                output.write(read)
        output.close()
        bam.close()

def valid_read(read, chrom, DBs, sizes, sizes2):
    start = read.reference_start
    end = read.reference_end
    span = end - start
    if span < sizes[0]:
        return False
    index = numpy.searchsorted(sizes2, span, side='left')
    n = span - sizes[index] + 1
    end = start + n
    if start >= end:
        return False
    try:
        scores = DBs[sizes[index]].values(chrom, start, end)
    except:
        return False
    if 1 in scores:
        return True
    else:
        return False


if __name__ == "__main__":
    main()

import numpy as np
import pandas as pd
from pyfaidx import Fasta

genomeFasta = "/zata/zippy/andrewsg/genome/T2T/chm13v1.fasta"

genome = Fasta(genomeFasta, 
               as_raw=True,
               sequence_always_upper=True)

with open("/zata/zippy/andrewsg/genome/T2T/chm13v1.chrom.sizes", "w") as f:
    for chrom in list(genome.keys()):
        print(chrom, len(genome[chrom][:]), file=f, sep="\t")
    



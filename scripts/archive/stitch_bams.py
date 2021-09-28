# Stitch two single end BAM files into a paired end bam file. 

## Setting up imports.
import pysam
import os
import shutil
import tempfile
import argparse
import sys

tmpdir=tempfile.mkdtemp(dir=".")
parser = argparse.ArgumentParser(description="Stitches two single end BAM files into a paired end bam file.")
parser.add_argument("-1", type=str, dest="bam1", help="path to BAM file containing the first read of the pair", required=True)
parser.add_argument("-2", type=str, dest="bam2", help="path to BAM file containing the second read of the pair", required=True)
parser.add_argument("-o", type=str, dest="output", help="path to the output BAM file", required=True)
args = parser.parse_args()

####################################################
################## RE-INTEGRATION ##################
####################################################
i=0
# Running through and re-integerating split and unsplit reads into a single file.
with pysam.Samfile(args.bam1, "rb") as sin1, \
        pysam.Samfile(args.bam2, "rb") as sin2, \
        pysam.Samfile(args.output, "wb", template=sin1) as sout:

    for read1 in sin1:
        try:
            read2=next(sin2)
        except StopIteration:
            sys.stdout.write("Ran out of end2 reads before end1 finished.")

        if read1.qname == read2.qname and read1.mapping_quality >= 30 and read2.mapping_quality >= 30:
            # We need to set appropriate flags for first/second read of a pair,
            # and for primary/secondary alignments. 5' reads are marked as primaries.
            read1.is_paired = read2.is_paired = True
            read1.is_read1 = True
            read2.is_read2 = True
            sout.write(read1)
            sout.write(read2)

        if i % 1000000 == 0:
            sys.stdout.write(f"Merging: {i} pairs merged\r")
        i += 1

sys.stdout.write(f"Sorting output file...")
bsorted=os.path.join(tmpdir, "sorted.bam")
pysam.sort("-o", bsorted, "-n", args.output)
shutil.move(bsorted, args.output)


# Mopping up.
shutil.rmtree(tmpdir)


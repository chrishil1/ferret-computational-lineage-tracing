# This file identifies unassigned reads. 

# To generate the input to this, I ran the following commands

# samtools sort -o barcode_headSorted_anno.bam barcode_headAligned_anno.bam

# samtools view \
#    -h barcode_headSorted_anno.bam | \
#    sed 's/Ferret-ASM1176430v1-1_//' | \
#    samtools view -bS - > reformatted_bam.bam

# featureCounts \
#    --extraAttributes 'gene_name,gene_biotype' \
#    -s 0 -p -B \
#    -t gene \
#    -a ~/Analysis/ref/annotation/ferret/ASM1176430v1.1/GCF_011764305.1_ASM1176430v1.1_genomic.gtf \
#    -R BAM \
#    -o feature_counts_test.txt \
#    reformatted_bam.bam

import pysam
from collections import defaultdict
import argparse
import os
import sys


def main():
    options = get_options()

    bin_dict = defaultdict(int)

    bin_size = 1000
    max_bin = 17080000
    all_bins = list(range(0, max_bin, bin_size))
    samfile = pysam.AlignmentFile(options.bam_file, "rb")
    for bin_position in all_bins:
        for read in samfile.fetch(options.chromosome, bin_position, bin_position + bin_size):
            all_tags = read.get_tags()
            for tag in all_tags:
                if tag[0] == "XS":
                    # Possible tags:
                    # Unassigned_MultiMapping
                    # Unasigned_NoFeatures
                    # Assigned
                    # Unassigned Ambiguity
                    if tag[1] == "Unassigned_NoFeatures":
                        bin_dict[bin_position] += 1
                        
    #print(bin_dict)
    sorted_dict = dict(sorted(bin_dict.items(), key=lambda item: item[1]))

    total_values = sum(bin_dict.values())

    with open(options.out_file, "w") as output:
        for i in sorted_dict:
            percent_reads = round(sorted_dict[i] / total_values * 100, 2)
            output.write("Position " + str(i) + ": " + str(sorted_dict[i]) + " reads, " + str(percent_reads) + " percent.\n")
    

def get_options():
    # House keeping to read in arguments from the command line
    parser = argparse.ArgumentParser()

    parser.add_argument("-b", "--bam_file", dest = "bam_file",
        help = "The bam file created by featurecounts using the -R flag.",
        action = "store",
        metavar = "\b")

    parser.add_argument("-o", "--output_file", dest = "out_file",
        help = "What file to write the results to.",
        action = "store",
        metavar = "\b")

    parser.add_argument("-c", "--chromosome", dest = "chromosome",
        help = "The chromosome to test.",
        action = "store",
        metavar = "\b")

    parser.add_argument("-v", "--verbose",
        help = "Make the output more verbose",
        action = "store_true")

    args = parser.parse_args()

    if not os.path.exists(args.bam_file):
        sys.exit("Fastq directory does not exist!\n" +
            "You can add your own using the -f argument. See --help for more information")

    return(args)



if __name__ == "__main__":
    main()

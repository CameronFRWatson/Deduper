#!/usr/bin/env python

#--------------------------------------------------------------------------------------------------------------
# Reference-based PCR duplicate removal tool
# Author: Cameron Watson
# Last Updated: Nov 24, 2020
#--------------------------------------------------------------------------------------------------------------

import argparse
import re

#--------------------------------------------------------------------------------------------------------------
# USER INPUT
#--------------------------------------------------------------------------------------------------------------

def get_args():
    parser = argparse.ArgumentParser(description = "A reference-based PCR duplicate removal tool. Takes a sorted sam file \
        containing uniquely mapped reads. Retains only the first read if duplicate reads are found based on UMI/randomer, \
            chromosome, and position.",
    add_help = True)

    parser.add_argument("-f", "--file", 
    help = "A SAM file that has been sorted by RNAME then POS to be deduplicated. \
        UMIs or randomers should be at end of QNAME", required = True)

    parser.add_argument("-p", "--paired", action = 'store_const', const = True,
    help = "Specify if input file contains paired-end reads. default is single-end", required = False)

    parser.add_argument("-u", "--umi", 
    help = "File containing UMIs, unset if randomers", required = False)

    parser.add_argument("--keepDupes", action = 'store_const', const = True,
    help = "Specify to write duplicates out to separate sam file, otherwise duplicates are just removed",
    required = False)

    parser.add_argument("--correctUMI", action = 'store_const', const = True,
    help = "Specify to error-correct UMIs that are one nucleotide away from a known UMI; otherwise, \
        reads with unknown UMIs are discarded. When specified, reads with UMIs that are more than \
            one nucleotide away from a known UMI are discarded. Corrected UMIs will appear in the output \
                sam file with an asterisc next to them. Randomers are not error-corrected, and specifying \
                    correctUMI without an umi file will result in error", required = False)
    
    return parser.parse_args()

args = get_args()

#--------------------------------------------------------------------------------------------------------------
# CLASSES, METHODS, AND FUNCTIONS
#--------------------------------------------------------------------------------------------------------------

class SingleEndRead:

    def __init__(self,read):
        # the whole read for output
        self.read = read
        # isolating PCR duplicate identifiers
        split_read = re.split("\t", read)
        self.umi = re.split(":", split_read[0])[-1] # this is a randomer if --umi not specified
        self.flag = int(split_read[1])
        self.chrom = split_read[2]
        self.pos = int(split_read[3])
        self.cig = split_read[5]

    def adjustFwd(self):
        '''
        Adjusts the left-most starting position for forward reads by subracting the amount of soft-clipping 
        that occurs on the left side of the read based on the CIGAR string
        '''
        # identify leftmost softclipping, then subtract value from pos
        split_cig = re.split('S', self.cig, maxsplit = 1)
        if split_cig[0].isnumeric():
            self.pos -= int(split_cig[0])
        return self

    def adjustRv(self):
        '''
        Adjusts the left-most starting position for reverse reads by adding the number of matches/mismatches, 
        insertions relative to the reference, and soft-clipping on the right side of the read based on 
        the CIGAR string
        '''
        # parse cigar string into individual components
        letters = []
        for char in self.cig:
            if char.isnumeric() == False:
                letters.append(char)
        values = re.split("[a-zA-Z]", self.cig)
        # a list of tuples containing each letter and its value in order
        value_pairs = zip(letters, values)
        pair_ctr = 0
        for pair in value_pairs:
            # skipping soft-clipping on the left-side
            if pair_ctr != 0:
                if pair[0] == "S":
                    # adding for soft-clipping on the right side
                    self.pos += int(pair[1])
            # adjusting for the length of the query sequence
            elif pair[0] == "M" or pair[0] == "I":
                self.pos += int(pair[1])
        return self

    def gen_dictKey(self):
        '''
        Generates a convenient dictionary key for a sam line based on it's components
        '''
        self.dkey = self.strand + self.umi + self.chrom + str(self.pos)
        return self

# UMI correction function
# -------------------------

def errCorrect_umi(x, udict):
    '''
    Calculates the Hamming Distance for a given incorrect UMI (x) against all UMIs in the specified UMI 
    dictionary (udict). If the Hamming Distance is equal to 1 for a single dictionary UMI, x is corrected
    to that UMI. If multiple known UMIs have a Hamming Dist equal to 1, x is not corrected and is discarded.
    Must specify --correctUMI for this function to be used.
    '''
    hd_list = []
    for key in udict:
        hamm_dist = 0
        for i in range(0, len(key)):
            if x[i] != key[i]:
                hamm_dist += 1 # calculating the hamming distance between unknown umi and current dict umi
        if hamm_dist == 1:
            hd_list.append(key) # storing the known UMIs that were close to the unknown
    if len(hd_list) == 1:
        return hd_list[0] # if the unknown is only close to a single known, return the known 
    else:
        return x

#--------------------------------------------------------------------------------------------------------------
# MAIN
#--------------------------------------------------------------------------------------------------------------

# Initial input checks
# -----------------------------

if args.paired:
    raise Exception("Paried end option not available yet, coming soon.")

# making sure that an UMI file is specified if error-correcting UMIs
if args.correctUMI:
    if not args.umi:
        raise Exception("Cannot specify correctUMI without UMI file. see --umi")

# creation of an UMI dictionary, skipped if randomers
if args.umi:

    umi_dict = {} # keys: unique UMIs, vals: times seen

    with open(args.umi, "r") as fh:
        for line in fh:
            umi = line.strip("\n")
            if umi in umi_dict:
                raise ValueError("Repeated UMIs present in UMI file: Remove repeats")
            else:
                umi_dict[umi] = 0

# Opening files
# --------------

# input sam file
sam_fh = open(args.file, "r")

# output file: deduplicated sam file
parse_path = re.split("/", args.file)[-1]
output_prefix = re.split(".sam", parse_path)[0]
output_name = output_prefix + "_deduped.sam"
output_fh = open(output_name, "w")

# optional duplicates file if specified
if args.keepDupes:
    dupfile_name = output_prefix + "_PCRduplicates.sam"
    dupe_fh = open(dupfile_name, "w")

# Initializing counters
# ------------------------------------------

dupli_counter = 0
unknown_umi_count = 0
read_counter = 0
corrected_umis = 0
umi_discards = 0
previous_chrom = 0 # just for the first instance of comparing chrom number

# Checking for and filtering PCR duplicates
# ------------------------------------------

for line in sam_fh:

    if re.match("^@", line):
        # write out the line to retain metadata
        output_fh.write(line)
    else:
        # make the line an instance of class SingleEndRead
        current_read = SingleEndRead(line)
        read_counter += 1
        # creates/clears a dictionary onced finished with a chromosome
        if current_read.chrom != previous_chrom:
            dupli_dict = {} # key: identifying info for reads, val: num of PCR duplicates with that info

        # check to see if line has a valid UMI, if --umi file given
        if args.umi:
            if current_read.umi not in umi_dict:
                unknown_umi_count +=1
                if args.correctUMI:
                    unknown_umi = current_read.umi
                    # updating with potentially corrected umi
                    current_read.umi = errCorrect_umi(unknown_umi, umi_dict)
                    # also updating the full read so it gets outputted with the corrected UMI
                    # UMI will have star next to it in QNAME if corrected from an unknown UMI
                    if current_read.umi != unknown_umi:
                        corrected_umis += 1
                        correction_flag = current_read.umi + "*"
                        current_read.read = re.sub(unknown_umi, correction_flag, current_read.read)
                else:
                    umi_discards += 1
                    continue # to the next line in the file, disregarding line with bad UMI
            if current_read.umi in umi_dict: 
                umi_dict[current_read.umi] += 1
            else:
                umi_discards += 1
                continue # to the next line in the file, disregarding line with bad UMI

        # check for strandedness from the bitwise flag and adjust pos accordingly
        if (current_read.flag & 16) == 16:
            current_read.strand = "r"
            current_read = current_read.adjustRv()
            current_read = current_read.gen_dictKey()
        else:
            current_read.strand = "f"
            current_read = current_read.adjustFwd()
            current_read = current_read.gen_dictKey()
        
        # checking current line against duplicate dictionary
        if current_read.dkey in dupli_dict:
            dupli_counter += 1
            dupli_dict[current_read.dkey] += 1
            if args.keepDupes:
                dupe_fh.write(current_read.read)
        else:
            dupli_dict[current_read.dkey] = 0
            output_fh.write(current_read.read)
        
        # tracking chromosome number so the dictionary can be cleared between chroms
        previous_chrom = current_read.chrom

# close all files
#-----------------

sam_fh.close()
output_fh.close()

if args.keepDupes:
    dupe_fh.close()

# summary file
summ_name = output_prefix + "_summary.txt"
with open(summ_name, "w") as fh:
    fh.write("Total number of reads:" + str(read_counter) + "\n")
    fh.write("Number of PCR duplicates removed:" + str(dupli_counter) + "\n")
    fh.write("Proportion of reads that were PCR duplicates:" + str(round(dupli_counter/read_counter, 2)) + "\n")
    if args.umi:
        fh.write("Number of unknown UMIs:" + str(unknown_umi_count) + "\n")
        if args.correctUMI:
            fh.write("Number of corrected UMIs:" + str(corrected_umis) + "\n")
            fh.write("Number of discarded UMIs:" + str(umi_discards) + "\n")


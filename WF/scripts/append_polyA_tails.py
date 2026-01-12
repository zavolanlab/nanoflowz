#!/usr/bin/env python3

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import warnings
warnings.simplefilter('ignore')

import sys
from argparse import ArgumentParser, RawTextHelpFormatter
import os
import subprocess
import pandas as pd
import numpy as np
import HTSeq
import pysam


def main():
    """ Parse bam file and append polyA tail of the length specified in the pt tag of each alignment
    """

    __doc__ = "append polyA tails to alignments, based on the pt tag"

    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument("--input_bam_file",
                        dest="input_bam_file",
                        help="Path to the input bam file",
                        required=True,
                        metavar="FILE",)
    parser.add_argument("--output_bam_file",
                        dest="output_bam_file",
                        help="Path to write the modified bam file",
                        required=True,
                        metavar="FILE",)
    try:
        options = parser.parse_args()
    except(Exception):
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    
    input_bam_file_path = options.input_bam_file
    output_bam_file_path = options.output_bam_file
    
    print('[INFO] processing bam file: '+input_bam_file_path+'\n')
    
    almnt_file = pysam.AlignmentFile(input_bam_file_path, "rb")
    bam_writer = pysam.AlignmentFile(output_bam_file_path, "wb", header=almnt_file.header)
    
    min_A_content = 0.9
    
    # all .iv elements follow bed format
    j,k=0,0
    total_number_of_reads,number_of_nonmodified_reads = 0,0
    for almnt in almnt_file:
        j=j+1
        pt = almnt.get_tag('pt')
        if almnt.is_forward:
            read_seq = almnt.query_sequence
            read_qualstr = almnt.query_qualities
            cigar = almnt.cigar
        else:
            read_seq = almnt.get_forward_sequence()
            read_qualstr = almnt.query_qualities[::-1]
            cigar = almnt.cigar[::-1]
            
        # find presumbale position of the cleavage site in the read - how many nt from 3'end to trim
        n = len(read_seq)
        Acount = 0
        suffix_len = 0
        prev_char = ''
        nonA_stretch_len = 0
        for i in range(n - 1, -1, -1):
            new_char = read_seq[i:i+1]
            suffix_len +=1
            if new_char.upper() == "A":
                Acount += 1
            else:
                if prev_char=="A":
                    nonA_stretch_len = 1
                else:
                    nonA_stretch_len += 1
                Aprop = Acount/suffix_len
                if Aprop<min_A_content:
                    break
            prev_char = new_char
        trim_len = suffix_len-nonA_stretch_len
    
        # modify sequence, qualseq and cigar
        if trim_len<pt:
            # define new read sequence
            new_read_seq = read_seq[:(-trim_len if trim_len>0 else None)]+"A"*pt
    
            # define new quality string
            quality_val = int(np.round(np.mean(read_qualstr[(-trim_len if trim_len>0 else -(min(5,len(read_seq)))):]),0)) # take the average of read qualities at trimmed part as a proxy read quality for added sequence
            new_read_qualstr = read_qualstr[:(-trim_len if trim_len>0 else None)]
            new_read_qualstr.extend([quality_val]*pt)
    
            # add some basic quality checks
            if len(new_read_seq)!=len(new_read_qualstr):
                print("[WARNING] modified read seq and read quality string have different length: "+str(almnt.query_name)+"; read is skipped")
                continue          
            
            # define new cigar string
            l=0 # position in cigar, backwards
            i = 0
            success = False
            for elem in cigar[::-1]:
                i=i+1
                if elem[0] in [0,1,4,5]:
                    l = l+elem[1]            
                    if l>trim_len: # we reached the right element
                        last_cigar_elem = elem
                        if last_cigar_elem[0]==4: # if soft-clip then extend soft-clip
                            cigar = cigar[:-i] # cigar until previous element will be untouched
                            cigar.append((4,last_cigar_elem[1]-trim_len+pt))
                            success=True
                        elif last_cigar_elem[0]==0: # if match then append soft-clip
                            cigar = cigar[:(-i+1 if -i+1<0 else None)]
                            cigar.append((4,(l-last_cigar_elem[1])+pt-trim_len))
                            success=True
                        else:
                            # this is undefined situation
                            # skip such read
                            success = False
                        break
            if not success:
                print("[WARNING] CIGAR string of read couldn't be modified: "+str(almnt.query_name)+"; read is skipped")
                continue
            if len(new_read_seq)!=np.sum([elem[1] for elem in cigar if elem[0] in [0,1,4,5]]):
                print("[WARNING] modified CIGAR and read seq have different lengths: "+str(almnt.query_name)+"; read is skipped")
                continue            
                
            if almnt.is_forward:
                almnt.query_sequence = new_read_seq
                almnt.query_qualities = new_read_qualstr
                almnt.cigar = cigar
            else:
                new_read_seq = HTSeq.Sequence( str.encode(new_read_seq), "seq" )
                almnt.query_sequence = new_read_seq.get_reverse_complement().seq.decode("utf-8")
                almnt.query_qualities = new_read_qualstr[::-1]
                almnt.cigar = cigar[::-1]
                
            # add a tag specifying that the read was modifed
            almnt.set_tag('pa',1,'i')
            
        else:
            # not clear how many to trim
            # leave reads as they are then
            number_of_nonmodified_reads +=1
            almnt.set_tag('pa',0,'i')
        
        bam_writer.write( almnt )
        k=k+1
    
    total_number_of_processed_reads = k
    total_number_of_input_reads = j
    print("[INFO] total number of non-modified reads: "+str(number_of_nonmodified_reads))
    print("[INFO] total number of input reads: "+str(total_number_of_input_reads))
    print("[INFO] total number of successfully processed reads: "+str(total_number_of_processed_reads))
    bam_writer.close()
    
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!")
        sys.exit(1)
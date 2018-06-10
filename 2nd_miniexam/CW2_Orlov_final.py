# -*- coding: utf-8 -*-
"""
Created on Sat Jun  9 17:25:24 2018

@author: Orlov
"""
import argparse
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO

if __name__ == "__main__":
    
    parser=argparse.ArgumentParser(description='Script for finding contigs that align on genomes with given E-value')
    
    parser.add_argument('-i', '--input', help='path to input fasta file', metavar='string', type=str)
    
    parser.add_argument('-e', '--e_value', help='E-value threshold', metavar='float number', type=float)
    
    args=parser.parse_args()
    
    input_file=args.input
    e_value_threshold=args.e_value

    
    app_output_path="aligned.fasta"
    not_app_output_path="nonaligned.fasta"
    
    with open (app_output_path, 'w') as out_al_file, open(not_app_output_path, 'w') as out_nonal_file:
        for fasta in SeqIO.parse(input_file, "fasta"):
            query = NCBIWWW.qblast("blastn", "nt", fasta.seq, expect=e_value_threshold, format_type="XML")
            blast_result = NCBIXML.parse(query)
            for result in blast_result:
                if len(result.alignments)>0:
                    SeqIO.write(fasta, out_al_file, "fasta")
                elif len(result.alignments)==0:
                    SeqIO.write(fasta, out_nonal_file, "fasta")




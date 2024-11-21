#! /usr/bin/env python3

import argparse
import os
from Bio import SeqIO


def filter_fasta_by_locus_tags(fasta_file, locus_tags_file, output_file):
    # Read the locus tags into a set for quick lookup
    with open(locus_tags_file, 'r') as lt_file:
        locus_tags = set(lt_file.read().splitlines())

    # Filter the FASTA records and write matching ones to the output file
    with open(output_file, 'w') as out_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            if record.id in locus_tags:
                SeqIO.write(record, out_file, 'fasta')

def main():
    parser = argparse.ArgumentParser(description='Filter FASTA file by locus tags.')
    parser.add_argument('fasta_file', type=str, help='Path to the input FASTA file.')
    parser.add_argument('locus_tags_file', type=str, help='Path to the file containing the locus tags.')
    parser.add_argument('output_file', type=str, help='Path for the output FASTA file.')

    args = parser.parse_args()

    # Ensure the parent directories for the output file exist
    os.makedirs(os.path.dirname(args.output_file), exist_ok=True)

    # Call the function with arguments
    filter_fasta_by_locus_tags(args.fasta_file, args.locus_tags_file, args.output_file)

if __name__ == "__main__":
    main()

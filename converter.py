#!/usr/bin/env

from Bio import SeqIO
import argparse

def convert_gb_to_fasta(gb_file, fasta_file):
    with open(fasta_file, 'w') as output_handle:
        for seq_record in SeqIO.parse(gb_file, "genbank"):
            organism = seq_record.annotations.get("organism", "unknown_organism").replace(" ", "_")
            genebank_id = seq_record.id.split('.')[0]
            for feature in seq_record.features:
                if feature.type == "CDS":
                    gene_seq = feature.extract(seq_record.seq)
                    gene_id = feature.qualifiers.get("gene", feature.qualifiers.get("product", ["unknown_gene"]))[0]
                    output_handle.write(f">{organism}_{genebank_id}_{gene_id}\n{gene_seq}\n")
                if feature.type == "rRNA":
                    gene_seq = feature.extract(seq_record.seq)
                    gene_id = feature.qualifiers.get("gene", feature.qualifiers.get("product", ["unknown_gene"]))[0]
                    output_handle.write(f">{organism}_{genebank_id}_{gene_id}\n{gene_seq}\n")


def main():
    parser = argparse.ArgumentParser(description="Convert GenBank file to FASTA format.")
    parser.add_argument("-g", "--gb", help="Path to the input GenBank file.", required=True)
    parser.add_argument("-o", "--output", help="Path to the output FASTA file.", required=True)
    args = parser.parse_args()

    convert_gb_to_fasta(args.gb, args.output)


if __name__ == '__main__':
    main()
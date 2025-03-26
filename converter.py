from Bio import SeqIO
import argparse

def convert_gb_to_fasta(gb_file, fasta_file):
    with open(fasta_file, 'w') as output_handle:
        for seq_record in SeqIO.parse(gb_file, "genbank"):
            for feature in seq_record.features:
                if feature.type == "gene":
                    gene_seq = feature.extract(seq_record.seq)
                    gene_id = feature.qualifiers.get("gene", feature.qualifiers.get("locus_tag", ["unknown_gene"]))[0]
                    output_handle.write(f">{gene_id}\n{gene_seq}\n")

def main():
    parser = argparse.ArgumentParser(description="Convert GenBank file to FASTA format.")
    parser.add_argument("gb", help="Path to the input GenBank file.")
    parser.add_argument("output", help="Path to the output FASTA file.")
    args = parser.parse_args()

    convert_gb_to_fasta(args.gb_file, args.fasta_file)


if __name__ == '__main__':
    main()
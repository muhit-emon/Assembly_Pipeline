import argparse
import pandas as pd 

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-aF", "--diamond_alignment_file", type = str, required = True, help = "Path to DIAMOND alignment file of Prodigal genes against DeepARG-DB or mobileOG-DB(required)")
    parser.add_argument("-gT", "--gene_type", type = str, required = True, help = "Alignment of which gene (ARG or MGE) (required)")
    parser.add_argument("-out", "--output_fname", type = str, required = True, help = "User provided output file name (required)")
    args = parser.parse_args()

    alignment_file = args.diamond_alignment_file
    gene_type = args.gene_type
    output_fname = args.output_fname

    filtered_diamond_output_file = output_fname + "_aln_filtered_" + gene_type + ".tsv"

    diamond_output = pd.read_csv(alignment_file, sep='\t', names=['qtitle','stitle', 'pident', 'bitscore', 'evalue', 'alignLen'], header=None)
    #diamond_output = diamond_output[diamond_output.bitscore >= 50] # bit-score filtering
    diamond_output = diamond_output[diamond_output.alignLen >= 25] # alignment length filtering

    diamond_output.to_csv(filtered_diamond_output_file, sep='\t', index=False)
    

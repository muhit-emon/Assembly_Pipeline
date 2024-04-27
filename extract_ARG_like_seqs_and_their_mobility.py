from Bio import SeqIO
import argparse
import pandas as pd 

def protein_dictionary(protein_fasta):

    ID_maps_seq = {}

    for sequence in SeqIO.parse(protein_fasta, "fasta"):
        this_prot_description = str(sequence.description)

        # extract ID from prodigal output
        last_part_of_header = this_prot_description.split("#")[-1].strip()
        ID_part = last_part_of_header.split(";")[0].strip()
        ID = ID_part.split("=")[-1].strip()

        this_prot = str(sequence.seq)
        ID_maps_seq[ID] = this_prot

    return ID_maps_seq

def output_ARG_like_seqs(ARG_alignment_file, ID_maps_seq, output_file):

    diamond_output = pd.read_csv(ARG_alignment_file, sep='\t', names=['qtitle','stitle', 'pident', 'bitscore', 'evalue'], header=None)

    f = open(output_file, "a") # output_fname + "_ARGs.faa"

    for i in range(len(diamond_output)):

        prodigal_header = diamond_output.iloc[i]["qtitle"]
        # extract ID from prodigal output
        last_part_of_header = prodigal_header.split("#")[-1].strip()
        ID_part = last_part_of_header.split(";")[0].strip()
        ID = ID_part.split("=")[-1].strip()

        t = prodigal_header.split("#")
        gene_location_info = t[0].strip() + " # " + t[1].strip() + " # " + t[2].strip() + " # " + t[3].strip() # something like scaffold_3_5 # 3697 # 4401 # -1

        DeepARG_DB_header = diamond_output.iloc[i]["stitle"]

        seq = ID_maps_seq[ID]

        f.write(">" + gene_location_info+"|"+DeepARG_DB_header)
        f.write("\n")
        f.write(seq)
        f.write("\n")
    
    f.close()

    return


def output_mobility_info(ARG_alignment_file, mobileOG_alignment_file, output_file):

    contigs_of_MGEs = {}
    diamond_output_mobileOG = pd.read_csv(mobileOG_alignment_file, sep='\t', names=['qtitle','stitle', 'pident', 'bitscore', 'evalue'], header=None)

    for i in range(len(diamond_output_mobileOG)):

        prodigal_header = diamond_output_mobileOG.iloc[i]["qtitle"]
        # extract ID from prodigal output
        last_part_of_header = prodigal_header.split("#")[-1].strip()
        ID_part = last_part_of_header.split(";")[0].strip()
        ID = ID_part.split("=")[-1].strip()

        contig_ID = ID.split("_")[0]

        contigs_of_MGEs[contig_ID] = 1
    

    diamond_output_ARG = pd.read_csv(ARG_alignment_file, sep='\t', names=['qtitle','stitle', 'pident', 'bitscore', 'evalue'], header=None)

    f = open(output_file, "a") # output_fname + "_ARGs_and_mobility.tsv"
    f.write("ARG-Header"+"\t"+"Co-occurs_with_MGE/Is_potentially_mobile"+"\n")

    for i in range(len(diamond_output_ARG)):

        prodigal_header = diamond_output_ARG.iloc[i]["qtitle"]
        # extract ID from prodigal output
        last_part_of_header = prodigal_header.split("#")[-1].strip()
        ID_part = last_part_of_header.split(";")[0].strip()
        ID = ID_part.split("=")[-1].strip()

        contig_ID = ID.split("_")[0]
        if contig_ID in contigs_of_MGEs:
            is_potentially_mobile = "YES"
        else:
            is_potentially_mobile = "NO"
        
        t = prodigal_header.split("#")
        gene_location_info = t[0].strip() + " # " + t[1].strip() + " # " + t[2].strip() + " # " + t[3].strip() # something like scaffold_3_5 # 3697 # 4401 # -1

        DeepARG_DB_header = diamond_output_ARG.iloc[i]["stitle"]

        ARG_header = gene_location_info+"|"+DeepARG_DB_header

        f.write(ARG_header + "\t" + is_potentially_mobile + "\n")

    
    f.close()

    return

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-pG", "--predicted_genes", type = str, required = True, help = "Prodigal predicted genes faa file (required)")
    parser.add_argument("-ARG_aF", "--ARG_alignment_file", type = str, required = True, help = "Path to DIAMOND alignment file of Prodigal genes against DeepARG-DB(required)")
    parser.add_argument("-mobileOG_aF", "--mobileOG_alignment_file", type = str, required = True, help = "Path to DIAMOND alignment file of Prodigal genes against mobileOG-DB(required)")
    parser.add_argument("-out", "--output_fname", type = str, required = True, help = "User provided output file name (required)")
    args = parser.parse_args()

    predicted_genes = args.predicted_genes
    ARG_alignment_file = args.ARG_alignment_file
    mobileOG_alignment_file = args.mobileOG_alignment_file
    output_fname = args.output_fname

    ID_maps_seq = protein_dictionary(predicted_genes)

    out_file_name_ARGs = output_fname + "_ARGs.faa"
    output_ARG_like_seqs(ARG_alignment_file, ID_maps_seq, out_file_name_ARGs)

    out_file_name_ARGs_and_their_mobility = output_fname + "_ARGs_and_mobility.tsv"
    output_mobility_info(ARG_alignment_file, mobileOG_alignment_file, out_file_name_ARGs_and_their_mobility)
import argparse
import pandas as pd 
from Bio import SeqIO

def protein_dictionary(protein_fasta):

    ID_maps_seq = {}

    for sequence in SeqIO.parse(protein_fasta, "fasta"):

        this_prot_description = str(sequence.description)

        # extract ID (e.g. contig_36_7) from prodigal output
        this_prodigal_header = this_prot_description
        tmp = this_prodigal_header.split(" # ")

        this_gene_prodigal_ID = tmp[0]

        this_prot = str(sequence.seq)
        ID_maps_seq[this_gene_prodigal_ID] = this_prot

    return ID_maps_seq

def output_ARG_like_seqs(ARG_alignment_file, ID_maps_seq, output_file):

    diamond_output = pd.read_csv(ARG_alignment_file, sep='\t')

    f = open(output_file, "a") # output_fname + "_ARGs.faa"

    for i in range(len(diamond_output)):

        prodigal_header = diamond_output.iloc[i]["qtitle"]

        # extract ID (e.g. contig_36_7) from prodigal output
        tmp = prodigal_header.split(" # ")
        this_gene_prodigal_ID = tmp[0]

        DeepARG_DB_header = diamond_output.iloc[i]["stitle"]

        seq = ID_maps_seq[this_gene_prodigal_ID]

        f.write(">" + DeepARG_DB_header)
        f.write("\n")
        f.write(seq)
        f.write("\n")
    
    f.close()

    return

def find_ARG_and_MGE_on_contigs(alignment_file):

    '''
    output: a dictionary with contig_ID as key and value is a list such as [(ABE02098|MLS|mefA|mefA, 36, 254), (BAL43359.1|MLS|mphG|mphG, 125, 786)]
    Here, (gene_info, gene_coordinate_start, gene_coordinate_end)
    '''
    d = {}

    df = pd.read_csv(alignment_file, sep='\t') # alignment_file has already headers ["qtitle", "stitle", "pident", "bitscore", "evalue", "alignLen"]

    for i in range(len(df)):

        this_prodigal_header = df.iloc[i]["qtitle"]
        tmp = this_prodigal_header.split(" # ")

        this_gene_prodigal_ID = tmp[0]
        t = this_gene_prodigal_ID.split("_")
        this_contig_ID = t[0]+"_"+t[1]

        this_gene_coordinate_start = int(tmp[1])
        this_gene_coordinate_end = int(tmp[2])

        this_gene_information = df.iloc[i]["stitle"]

        if this_contig_ID not in d:
            d[this_contig_ID] = [(this_gene_information, this_gene_coordinate_start, this_gene_coordinate_end)]
        else:
            d[this_contig_ID].append((this_gene_information, this_gene_coordinate_start, this_gene_coordinate_end))

    return d

def read_k2_output_and_create_dictionary(k2_output_file):

    # Read the TSV file without header
    df = pd.read_csv(k2_output_file, sep='\t', names=['is_classified','seq_ID', 'taxa', 'len', 'LCA'], header=None)

    # Filter the dataframe if the contigs are "C"/classified
    filtered_df = df[df['is_classified'] == 'C']

    # Create a dictionary using the contig_ID as keys and its taxa as value
    d = dict(zip(filtered_df['seq_ID'], filtered_df['taxa']))
    return d

def mobileOG_DB_MGE_class_dictionary(mobileOG_DB_class_tsv_file):

    df = pd.read_csv(mobileOG_DB_class_tsv_file, sep="\t")
    mobileOG_entry_name_maps_MGE_class = dict(zip(df['mobileOG_Entry_Name'], df['MGE_Class'].str.lower()))
    return mobileOG_entry_name_maps_MGE_class

def output_ARGs_and_their_cooccurrence_with_MGEs_and_pathogens(contig_ID_maps_taxa, contig_ID_maps_ARGs, contig_ID_maps_MGEs, mobileOG_ID_maps_MGE_class, output_file):

    f = open(output_file, "a") # output_fname + "_ARGs_and_mobility_and_co-occurrence_with_pathogens.tsv"
    f.write("Accession|Drug|ARG"+"\t"+"Does_co-occur_with_MGE?"+"\t"+"Pathogen_co-occurrence"+"\n")

    for contig in contig_ID_maps_ARGs:
        this_contig_ARG_list = contig_ID_maps_ARGs[contig] # [(ABE02098|MLS|mefA|mefA, 36, 254), (BAL43359.1|MLS|mphG|mphG, 125, 786)]. Here, (gene_info, gene_coordinate_start, gene_coordinate_end)

        if contig in contig_ID_maps_MGEs:
            this_contig_MGE_list = contig_ID_maps_MGEs[contig]
        else:
            this_contig_MGE_list = []

        if contig in contig_ID_maps_taxa:
            pathogen_flag = True
        else:
            pathogen_flag = False

        for (ARG, ARG_gc_start, ARG_gc_end) in this_contig_ARG_list:
            max_distance_bw_ARG_and_MGE_to_get_mobile_ARG = 10000 # maximum distance between ARG and MGE is set to 10kb to call an ARG as mobile ARG
            '''
            for only "plasmid" tag in mobileOG_DB class information, we will just check whether an ARG and that plasmid co-occurs on the same contig
            for any other class tag in mobileOG_DB class information, we will check if the ARG/MGE distance is <= 10kbp
            '''
            MGE_flag = False
            for (MGE, MGE_gc_start, MGE_gc_end) in this_contig_MGE_list:
                this_MGE_ID = MGE.split("|")[0]
                if this_MGE_ID in mobileOG_ID_maps_MGE_class:
                    this_MGE_class = mobileOG_ID_maps_MGE_class[this_MGE_ID]
                else:
                    this_MGE_class = "none" # not available in mobileOG MGE class metadata
                    
                if this_MGE_class == "plasmid":
                    MGE_flag = True
                else:
                    distances = [abs(ARG_gc_start-MGE_gc_start), abs(ARG_gc_start-MGE_gc_end), abs(ARG_gc_end-MGE_gc_start), abs(ARG_gc_end-MGE_gc_end)]
                    distance = min(distances)  # Smallest distance between any endpoints
                    if distance <= max_distance_bw_ARG_and_MGE_to_get_mobile_ARG:
                        MGE_flag = True

            this_ARG_accession = ARG.split("|")[0]
            this_ARG_drug = ARG.split("|")[1]
            this_ARG_gene_name = ARG.split("|")[2]
            this_ARG_accession_drug_and_gene_name = this_ARG_accession + "|" + this_ARG_drug + "|" + this_ARG_gene_name

            if MGE_flag is True and pathogen_flag is True:
                f.write(this_ARG_accession_drug_and_gene_name+"\t"+"YES"+"\t"+contig_ID_maps_taxa[contig]+"\n")
            elif MGE_flag is True and pathogen_flag is False:
                f.write(this_ARG_accession_drug_and_gene_name+"\t"+"YES"+"\t"+"N/A"+"\n")
            elif MGE_flag is False and pathogen_flag is True:
                f.write(this_ARG_accession_drug_and_gene_name+"\t"+"NO"+"\t"+contig_ID_maps_taxa[contig]+"\n")
            elif MGE_flag is False and pathogen_flag is False:
                f.write(this_ARG_accession_drug_and_gene_name+"\t"+"NO"+"\t"+"N/A"+"\n")
    f.close()
    return

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-pP", "--predicted_proteins", type = str, required = True, help = "Prodigal predicted proteins (translated) faa file (required)")
    parser.add_argument("-ARG_aF", "--ARG_alignment_file", type = str, required = True, help = "Path to DIAMOND alignment file of Prodigal genes against DeepARG-DB(required)")
    parser.add_argument("-mobileOG_aF", "--mobileOG_alignment_file", type = str, required = True, help = "Path to DIAMOND alignment file of Prodigal genes against mobileOG-DB(required)")
    parser.add_argument("-mobileOG_class_info", "--mobileOG_class_info", type = str, required = True, help = "Path to mobileOG-DB MGE class information (required)")
    parser.add_argument("-k2_output", "--k2_output", type = str, required = True, help = "Kraken2 contig taxonomy classification file (required)")
    parser.add_argument("-out", "--output_fname", type = str, required = True, help = "User provided output file name (required)")
    args = parser.parse_args()

    
    predicted_proteins = args.predicted_proteins
    ARG_alignment_file = args.ARG_alignment_file
    mobileOG_alignment_file = args.mobileOG_alignment_file
    mobileOG_class_info_file = args.mobileOG_class_info
    kraken2_output_file = args.k2_output
    output_fname = args.output_fname
    
    mobileOG_entry_name_maps_MGE_class = mobileOG_DB_MGE_class_dictionary(mobileOG_class_info_file)
    contig_ID_maps_taxa = read_k2_output_and_create_dictionary(k2_output_file=kraken2_output_file)
    contig_ID_maps_ARGs = find_ARG_and_MGE_on_contigs(alignment_file=ARG_alignment_file)
    contig_ID_maps_MGEs = find_ARG_and_MGE_on_contigs(alignment_file=mobileOG_alignment_file)

    detected_ARGs_and_co_occurrence_with_MGEs_and_pathogens_file = output_fname + "_ARGs_and_mobility_and_co-occurrence_with_pathogens.tsv"
    detected_ARGs_output_file = output_fname + "_ARGs.faa"

    output_ARGs_and_their_cooccurrence_with_MGEs_and_pathogens(contig_ID_maps_taxa, contig_ID_maps_ARGs, contig_ID_maps_MGEs, mobileOG_entry_name_maps_MGE_class, detected_ARGs_and_co_occurrence_with_MGEs_and_pathogens_file)

    ID_maps_seq = protein_dictionary(predicted_proteins)
    output_ARG_like_seqs(ARG_alignment_file, ID_maps_seq, detected_ARGs_output_file)
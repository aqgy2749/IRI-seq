import os
import numpy as np
import pandas as pd
import gzip
import sys
import pickle
import csv
import shutil

def Gene_count_summary(input_folder, sample_ID_file, output_folder):

    sample_ID = list(pd.read_csv(sample_ID_file, header=None)[0])
    #input_files = [input_folder + "/" + sample_name + ".count.gz" for sample_name in sample_ID]
    if os.path.exists(output_folder):
        print("Output folder exists. Delete the current folder.")
        shutil.rmtree(output_folder)
        
    print("Make new output folder: ", output_folder)
    os.makedirs(output_folder)

    Gene_count_output = os.path.join(output_folder, "Gene_count.txt.gz")
    Exon_count_output = os.path.join(output_folder, "Exon_count.txt.gz")
    Intron_count_output = os.path.join(output_folder, "Intron_count.txt.gz")
    Exon_ID_count_output = os.path.join(output_folder, "Exon_count_with_exon_id.txt.gz")
    
    for sample_name in sample_ID:
        input_file = input_folder + "/" + sample_name + ".count.gz"
        print("***********************Process input file: ", input_file)
        df = pd.read_csv(input_file, sep="\t", compression="gzip", header=None,index_col=False, dtype={0: str},
                        names = ["chrom", "start", "end", "read_name", "mapped", "strand", "Exon_intron", "Single_multiple_gene", "Gene_id",
                                "Nearest_gene_5prime", "Nearest_gene_3prime", "Exon_id", "Mapping_strand"])
        df[["barcode", "UMI"]] = df['read_name'].str.split(',', n=3, expand=True).iloc[:, :2]
        df["cell_name"] = sample_name + "." + df['barcode'].astype(str)
        # filter reads that map to a single gene
        #print("Number of reads before filtering reads that mapped to a single gene: ", df.shape[0])
        #df = df[df["Gene_id"] == "Single_gene"]
        print("Percentage of reads that mapped to single gene: ", df[df["Single_multiple_gene"] == "Single_gene"].shape[0] / df.shape[0])
        print("Percentage of reads that mapped to multiple genes: ", df[df["Single_multiple_gene"] == "Multi_gene"].shape[0] / df.shape[0])
        print("Percentage of reads that mapped to no genes: ", df[df["Single_multiple_gene"] == "Not_aligned"].shape[0] / df.shape[0])

        print("Number of reads before filtering reads that mapped to multiple genes: ", df.shape[0])
        df = df[df["Single_multiple_gene"] != "Multi_gene"]
        print("Number of reads after filtering reads that mapped to multiple genes: ", df.shape[0])

        # output the gene count matrix
        print("Process and output the gene count matrix...")
        df_gene_count = df.groupby(["cell_name", "Nearest_gene_3prime"])["UMI"].nunique().reset_index()
        with gzip.open(Gene_count_output, "at") as f:
            df_gene_count.to_csv(f, sep="\t", index=False, header = None)
        # output the exon count matrix
        print("Process and output the exon count matrix...")
        df_gene_count = df[df["Exon_intron"] == "Exon"].groupby(["cell_name", "Nearest_gene_3prime"])["UMI"].nunique().reset_index()
        with gzip.open(Exon_count_output, "at") as f:
            df_gene_count.to_csv(f, sep="\t", index=False, header = None)
        # output the intron count matrix
        print("Process and output the intron count matrix...")
        df_gene_count = df[df["Exon_intron"] == "Intron"].groupby(["cell_name", "Nearest_gene_3prime"])["UMI"].nunique().reset_index()
        with gzip.open(Intron_count_output, "at") as f:
            df_gene_count.to_csv(f, sep="\t", index=False, header = None)
        # output the exon matrix by exon ID
        print("Process and output the exon count matrix by exon ID...")
        df_gene_count = df[(df["Exon_id"] != "No_unique_exon_match") & (df["Single_multiple_gene"] == "Single_gene") & (df["Exon_intron"] == "Exon")].groupby(["cell_name", "Exon_id"])["UMI"].nunique().reset_index()
        with gzip.open(Exon_ID_count_output, "at") as f:
            df_gene_count.to_csv(f, sep="\t", index=False, header = None)

    print("All input files are processed.")
'''
if __name__ == "__main__":
    input_folder = sys.argv[1]
    sample_ID = sys.argv[2]
    output_folder = sys.argv[3]
    Gene_count_summary(input_folder, sample_ID, output_folder)
'''
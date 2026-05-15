import os
import numpy as np
import pandas as pd
import gzip
import sys
import pickle
import csv
import shutil
import scanpy as sc
from scipy.sparse import csr_matrix

# create sparase matrix
def make_sparse_matrix(df_gene_count):
    # Create a dictionary that maps gene names to column indices
    gene_names = sorted(df_gene_count["Gene_id"].unique())
    gene_to_col = {gene_name: i for i, gene_name in enumerate(gene_names)}

    # Create a dictionary that maps cell names to row indices
    cell_names = sorted(df_gene_count["Cell_name"].unique())
    cell_to_row = {cell_name: i for i, cell_name in enumerate(cell_names)}

    # Create an array for the values of the matrix
    values = df_gene_count["UMI_count"].values

    # Create arrays for the row and column indices of the matrix
    rows = np.array([cell_to_row[cell_name] for cell_name in df_gene_count["Cell_name"]])
    cols = np.array([gene_to_col[gene_name] for gene_name in df_gene_count["Gene_id"]])

    # Create the sparse matrix
    sparse_matrix = csr_matrix((values, (rows, cols)), shape=(len(cell_names), len(gene_names)))
    df_cell = pd.DataFrame({"Cell_name" : cell_names})
    df_gene = pd.DataFrame({"Gene_id" : gene_names})
    
    return([sparse_matrix, df_cell, df_gene])

def generate_adata_from_gene_count(input_folder, df_gene_file, output_folder, UMI_limit):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    df_gene_count_file = os.path.join(input_folder, "Gene_count.txt.gz")
    df_exon_count_file = os.path.join(input_folder, "Exon_count.txt.gz")
    df_intron_count_file = os.path.join(input_folder, "Intron_count.txt.gz")
    df_Exon_ID_count_file = os.path.join(input_folder, "Exon_count_with_exon_id.txt.gz")

    df_gene_count = pd.read_csv(df_gene_count_file, sep="\t", compression="gzip", header=None, index_col=False, dtype={0: str},
                            names = ["Cell_name", "Gene_id", "UMI_count"])

    # how many count per cell do we have (> UMI_limit)
    df_cell_UMI_count = df_gene_count.groupby("Cell_name")["UMI_count"].sum().reset_index().sort_values("UMI_count", ascending = False)
    print("Total number of cell barcode detected: ", len(df_cell_UMI_count))
    df_cell_filtered = df_cell_UMI_count[df_cell_UMI_count["UMI_count"] > UMI_limit]
    print("Total number of cell barcode with UMI count > UMI limit: ", len(df_cell_filtered))
    print("Fraction of UMIs that are from real cells: ", np.sum(df_cell_filtered["UMI_count"]) / np.sum(df_cell_UMI_count["UMI_count"]))
    df_gene_count = df_gene_count[df_gene_count["Cell_name"].isin(df_cell_filtered["Cell_name"])]

    # group by cell name and aggregate UMI counts
    df_cell_grouped = df_gene_count.groupby("Cell_name").agg({"UMI_count": "sum"})

    # count Not_aligned UMIs for each cell
    df_not_aligned = df_gene_count[df_gene_count["Gene_id"] == "Not_aligned"].groupby("Cell_name").agg({"UMI_count": "sum"})

    # merge the two dataframes on cell name
    df_cell_grouped = df_cell_grouped.merge(df_not_aligned, on="Cell_name", how="left")

    # calculate the ratio of not_aligned UMIs to total UMIs
    df_cell_grouped["not_aligned_ratio"] = df_cell_grouped["UMI_count_y"] / df_cell_grouped["UMI_count_x"]

    # rename columns
    df_cell_grouped = df_cell_grouped.rename(columns={"UMI_count_x": "total_UMI_count", "UMI_count_y": "not_aligned_UMI_count"}).reset_index()
    df_gene_count = df_gene_count[df_gene_count["Gene_id"] != "Not_aligned"]

    df_result = make_sparse_matrix(df_gene_count)
    df_cell_full = df_result[1]
    df_gene_full = df_result[2]
    df_cell_full = df_cell_full.merge(df_cell_grouped, on="Cell_name", how="left")


    # Read the exonic count matrix
    df_gene_count = pd.read_csv(df_exon_count_file, sep="\t", compression="gzip", header=None, index_col=False, dtype={0: str},
                            names = ["Cell_name", "Gene_id", "UMI_count"])
    df_gene_count = df_gene_count[df_gene_count["Cell_name"].isin(df_cell_full["Cell_name"])]
    df_cell_grouped = df_gene_count.groupby("Cell_name").agg({"UMI_count": "sum"}).reset_index().rename(columns={"UMI_count" : "Exon_count"})
    df_cell_full = df_cell_full.merge(df_cell_grouped, on="Cell_name", how="left")
    df_cell_full["Exon_ratio"] = df_cell_full.Exon_count / df_cell_full.total_UMI_count
    df_cell_full["UMI_count"] = df_cell_full.total_UMI_count - df_cell_full.not_aligned_UMI_count
    df_result_exon = make_sparse_matrix(df_gene_count)

    # Read the intronic count matrix
    df_gene_count = pd.read_csv(df_intron_count_file, sep="\t", compression="gzip", header=None, index_col=False, dtype={0: str},
                            names = ["Cell_name", "Gene_id", "UMI_count"])
    df_gene_count = df_gene_count[df_gene_count["Cell_name"].isin(df_cell_full["Cell_name"])]
    df_result_intron = make_sparse_matrix(df_gene_count)

    # Annotate df_gene
    df_gene = pd.read_csv(df_gene_file, sep = ",", index_col=False)
    df_gene.columns = ["Gene_id", "Gene_type", "Gene", "Gene_name", "Index_ID"]
    print("Any gene id not in the gene reference list: ", np.sum([not x for x in list(df_gene_full.Gene_id.isin(df_gene.Gene_id))]))

    df_gene_full = df_gene_full.merge(df_gene[["Gene_id", "Gene_type", "Gene_name"]], on = "Gene_id", how = "left")


    # Generate adata for full gene count matrix
    adata_full = sc.AnnData(X = df_result[0])
    df_cell_full.index = df_cell_full.Cell_name
    df_cell_full = df_cell_full.drop(columns = ["Cell_name"])
    df_gene_full.index = df_gene_full.Gene_id
    df_gene_full = df_gene_full.drop(columns = ["Gene_id"])
    adata_full.obs = df_cell_full
    adata_full.var = df_gene_full

    # Create adata for exonic count
    adata_exon = sc.AnnData(X = df_result_exon[0])
    df_cell_exon = df_result_exon[1]
    df_gene_exon = df_result_exon[2]
    df_cell_exon.index = df_cell_exon.Cell_name
    df_gene_exon.index = df_gene_exon.Gene_id

    df_cell_exon = df_cell_exon.drop(columns = ["Cell_name"])
    df_gene_exon = df_gene_exon.drop(columns = ["Gene_id"])

    adata_exon.obs = df_cell_exon
    adata_exon.var = df_gene_exon

    # Create adata for intronic count matrix
    adata_intron = sc.AnnData(X = df_result_intron[0])
    df_cell_intron = df_result_intron[1]
    df_gene_intron = df_result_intron[2]
    df_cell_intron.index = df_cell_intron.Cell_name
    df_gene_intron.index = df_gene_intron.Gene_id

    df_cell_intron = df_cell_intron.drop(columns = ["Cell_name"])
    df_gene_intron = df_gene_intron.drop(columns = ["Gene_id"])

    adata_intron.obs = df_cell_intron
    adata_intron.var = df_gene_intron

    print("Writing gene count matrix, exon count matrix and intron count matrix.....")
    adata_full.write_h5ad(os.path.join(output_folder, "adata_full.h5ad"))
    adata_exon.write_h5ad(os.path.join(output_folder, "adata_exon.h5ad"))
    adata_intron.write_h5ad(os.path.join(output_folder, "adata_intron.h5ad"))
'''
if __name__ == "__main__":
    input_folder = sys.argv[1]
    df_gene_file = sys.argv[2]
    output_folder = sys.argv[3]
    UMI_limit = int(sys.argv[4])
    generate_adata_from_gene_count(input_folder, df_gene_file, output_folder, UMI_limit)
'''
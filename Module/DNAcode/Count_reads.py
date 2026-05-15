import subprocess
import csv
import os
import pandas as pd
import gzip
import scanpy as sc
from multiprocessing import Pool

def Fastq_count_reads(fastq_file):
    with gzip.open(fastq_file, 'rt') as f:  # 'rt' mode opens the file in text mode
        count = sum(1 for line in f)
    return count // 4  # Divide by 4 for the four lines per read

def read_csv_to_list(file_name):
    return list(pd.read_csv(file_name, header=None)[0])

def Fastq_count_reads_files(input_folder, sample_ID, core = 5):
    sample_list = read_csv_to_list(sample_ID)
    with Pool(processes=core) as pool:
        result = pool.map(Fastq_count_reads, [(input_folder + sample + ".R2.fastq.gz") for sample in sample_list])
    return result

def SAM_count_mapped_reads(samfile_path):
    with open(samfile_path, 'r') as samfile:
        mapped_reads = 0
        for line in samfile:
            if line[0] != '@': # this line is not a header
                fields = line.strip().split('\t')
                flag = int(fields[1])
                if flag != 4:  # if the flag is not 4, the read is mapped
                    mapped_reads += 1
    return mapped_reads

def SAM_count_mapped_reads_files(input_folder, sample_ID, core = 5):
    sample_list = read_csv_to_list(sample_ID)
    with Pool(processes=core) as pool:
        result = pool.map(SAM_count_mapped_reads, [(input_folder + "/" + sample + "Aligned.out.sam") for sample in sample_list])
    return result

def Count_Align_STAR(input_file):
    df_star = pd.read_csv(input_file, header= None, sep = "\t")
    result = [int(i) for i in (list(df_star.iloc[[4, 7, 22], 1]))]
    return result
def Count_Align_STAR_files(input_folder, sample_ID, core = 5):
    sample_list = read_csv_to_list(sample_ID)
    with Pool(processes=core) as pool:
        result = pool.map(Count_Align_STAR, [(input_folder + "/" + sample + "Log.final.out") for sample in sample_list])
    result = pd.DataFrame(result, index = sample_list)
    return result


def SAM_count_reads(sam_file):
    with open(sam_file, 'r') as f:  
        count = sum(1 for line in f if not line.startswith('@'))
    return count

def SAM_count_reads_files(input_folder, sample_ID, core = 5):
    sample_list = read_csv_to_list(sample_ID)
    with Pool(processes=core) as pool:
        result = pool.map(SAM_count_reads, [(input_folder + "/" + sample + ".sam") for sample in sample_list])
    return result

def count_mapped_reads(gene_count_file):
    df_tmp = pd.read_csv(gene_count_file, sep = ",", header = None, index_col=False, names = ["Type", "Sample_name", "Number"])
    result = df_tmp.iloc[0:8, 2].sum()
    return(result)
def count_mapped_reads_files(input_folder, sample_ID, core = 5):
    sample_list = read_csv_to_list(sample_ID)
    with Pool(processes=core) as pool:
        result = pool.map(count_mapped_reads, [(input_folder + "/" + sample + ".report") for sample in sample_list])
    return result

def Count_cell_reads(Adata_folder):
    adata = sc.read_h5ad(Adata_folder + "/adata_full.h5ad")
    df_tmp = adata.obs
    df_tmp["Exp"] = [i.split(".")[0] for i in list(df_tmp.index)]
    df_grouped = df_tmp.groupby("Exp")["UMI_count"].sum().reset_index()
    return list(df_grouped["UMI_count"])
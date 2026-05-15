import subprocess
import csv
import os
from multiprocessing import Pool
from functools import partial

def Fastq_trim_sample(sample, input_folder, output_folder, cutadapt, adapter_seq):
    input_command = f"{cutadapt} -a {adapter_seq} --trim-n --minimum-length 20 -o {output_folder}/{sample}.R2.fastq.gz {input_folder}/{sample}.R2.fastq.gz"
    print(input_command)
    result = subprocess.check_output(input_command, shell=True, text=True)
    print(result)

def Fastq_trim_files(input_folder, sample_ID, output_folder, core, adapter_seq = "AAAAAAAA", cutadapt_path = "/ru-auth/local/home/jcao/anaconda3_new/envs/cutadaptenv/bin/cutadapt"):
    print("input folder:", input_folder)
    print("sample ID:", sample_ID)
    print("output folder:", output_folder)
    print("core number:", core)
    print("Adapter sequence:", adapter_seq)
    
    sample_file = open(sample_ID)
    sample_list = []
    for line in sample_file:
        sample = line.strip()
        sample_list.append(sample)
    sample_file.close()

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Update this part to pass the extra arguments to function_sample
    p = Pool(processes=int(core))
    func = partial(Fastq_trim_sample, input_folder = input_folder, output_folder=output_folder,
                  cutadapt = cutadapt_path, adapter_seq = adapter_seq)
    print("\n***************",sample_list)
    result = p.map(func, sample_list)
    p.close()
    p.join()

    print("All samples are processed.")

import subprocess
import csv
import os
from multiprocessing import Pool

def function_sample(sample, input_folder, output_folder, samtools_path="/rugpfs/fs0/cao_lab/scratch/asziraki/anaconda3/envs/original_pipeline/bin/samtools"):
    input_command = f"{samtools_path} view -bh -q 30 -F 4 {input_folder}/{sample}Aligned.out.sam|{samtools_path} sort -@ 10 -|{samtools_path} view -h ->{output_folder}/{sample}.sam"
    print(input_command)
    result = subprocess.check_output(input_command, shell=True, text=True)
    print(result)

    
def Sam_filter_files(input_folder, sample_ID, output_folder, core, samtools_path = "/rugpfs/fs0/cao_lab/scratch/asziraki/anaconda3/envs/original_pipeline/bin/samtools"):
    print("input folder:", input_folder)
    print("sample ID:", sample_ID)
    print("output folder:", output_folder)
    print("core number:", core)

    with open(sample_ID) as f:
        reader = csv.reader(f, delimiter='\t')
        sample_name = [row[0] for row in reader]

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Update this part to pass the extra arguments to function_sample
    with Pool(processes=core) as pool:
        pool.starmap(function_sample, [(name, input_folder, output_folder, samtools_path) for name in sample_name])

    print("All samples are processed.")

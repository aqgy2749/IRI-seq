import subprocess
import sys
from multiprocessing import Pool
from functools import partial
import pandas as pd
import gzip

def remove_duplicates(sample, input_folder, output_folder):
    # Read the zipped text file into a pandas DataFrame
    input_file = input_folder + "/" + sample + ".spatial.txt.gz"
    output_file = output_folder + "/" + sample + ".spatial.csv.gz"
    df = pd.read_csv(input_file, compression='gzip', sep='\t')

    # Drop duplicate rows based on all columns
    df_unique = df.drop_duplicates()

    # Write the unique rows to a new zipped text file
    df_unique.to_csv(output_file, sep='\t', compression='gzip', index=False)
    
    print("sample name: %s, total line: %f, filtered line: %f, filter rate: %f"
          % (sample, len(df), len(df_unique), float(len(df_unique)) / float(len(df))))

def remove_duplicates_files(input_folder, sampleID, output_folder, core):
    init_message = '''
    --------------------------start attaching UMI-----------------------------
    input folder: %s
    sample ID: %s
    output_folder: %s
    ___________________________________________________________________________
    ''' % (input_folder, sampleID, output_folder)

    print(init_message)


    # for each sample in the sample list, use the read1 file, read2 file, output file
    # and barcode_list to run UMI_attach_read2_barcode_list
    sample_file = open(sampleID)
    sample_list = []
    for line in sample_file:
        sample = line.strip()
        sample_list.append(sample)
    sample_file.close()

    # parallel for the functions
    p = Pool(processes=int(core))
    func = partial(remove_duplicates, input_folder=input_folder, output_folder=output_folder)
    result = p.map(func, sample_list)
    p.close()
    p.join()

    # print the completion message
    com_message = '''~~~~~~~~~~~~~~~Duplicates removal done~~~~~~~~~~~~~~~~~~'''
    print(com_message)

if __name__ == "__main__":
    input_folder = sys.argv[1]
    sampleID = sys.argv[2]
    output_folder = sys.argv[3]
    core = sys.argv[4]
    remove_duplicates_files(input_folder, sampleID, output_folder, core)

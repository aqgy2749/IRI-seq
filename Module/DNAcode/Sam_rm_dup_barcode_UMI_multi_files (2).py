
import sys
import pandas as pd
from collections import Counter
import subprocess
from multiprocessing import Pool
from functools import partial
from Levenshtein import distance
sys.setrecursionlimit(1000000)


def rm_dup_samfile(sample, input_folder, output_folder):
    samfile = input_folder + "/" + sample + ".sam"
    output_file = output_folder + "/" + sample + ".sam"
    
    f1 = open(samfile)
    f2 = open(output_file, 'w')

    pre_barcode = set()
    pre_line = []
    unique_id = []
    pre_chrom = 0
    pre_site = 0
    total_read_count = 0
    unique_read_count = 0


    for line in f1:
        
        if (line[0] == '@'):
            f2.write(line)
        else:
            total_read_count += 1
            fields = line.split('\t')
            name = fields[0].split(',')
            barcode_UMI = name[0] + name[1]
            chrom_num = fields[2]
            start_site = (fields[3])

            if ((start_site == (pre_site)) and (chrom_num == pre_chrom)):
                pre_chrom = chrom_num
                pre_site = start_site
                if barcode_UMI not in pre_barcode:
                    f2.write(line)
                    unique_read_count += 1
                    pre_barcode.add(barcode_UMI)

                else:
                    continue
            else:
                f2.write(line)
                unique_read_count += 1
                pre_chrom = chrom_num
                pre_site = start_site
                pre_barcode = set([barcode_UMI])
    
    f1.close()
    f2.close()
    
    end_message = '''
    ---------------------------------------------------------------------------
    sample ID: {0}
    Total input read number: {1}
    Total unique read number: {2}
    Duplication rate: {3}
    ___________________________________________________________________________
    '''.format(sample, str(total_read_count), str(unique_read_count), 1 - float(unique_read_count) / float(total_read_count))
    
    print(end_message)


# this function accept an input folder and a output folder and then generate the output file with the index
def rm_dup_files(input_folder, sampleID, output_folder, core):
    
    init_message = '''
    --------------------------start removing duplicates-----------------------------
    input folder: {0}
    sample ID: {1}
    output_folder: {2}
    ___________________________________________________________________________
    '''.format(input_folder, sampleID, output_folder)
    
    print(init_message)
    
    
    #for each sample in the sample list, remove duplicates for the input file in the input folder and 
    # output the sam file into a new folder
    
    sample_file = open(sampleID)
    sample_list = []
    for line in sample_file:
        sample = line.strip()
        sample_list.append(sample)
    sample_file.close()

    # parallele for the functions
    p = Pool(processes = int(core))

    func = partial(rm_dup_samfile, input_folder = input_folder, output_folder = output_folder)
    result = p.map(func, sample_list)
    p.close()
    p.join()
    
    #print the completion message
    com_message = '''~~~~~~~~~~~~~~~Duplicate removal done~~~~~~~~~~~~~~~~~~'''
    print(com_message)


if __name__ == "__main__":
    input_folder = sys.argv[1]
    sampleID = sys.argv[2]
    output_folder = sys.argv[3]
    core=sys.argv[4]
    rm_dup_files(input_folder, sampleID, output_folder, core)

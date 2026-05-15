import subprocess
import sys
from Levenshtein import distance
import gzip
from multiprocessing import Pool
from functools import partial
import pickle

def extract_spatial_barcode(sample, input_folder, output_folder, list_barcode_1, list_barcode_2, list_barcode_3, list_barcode_4, mismatch_rate):
    # open the read1, read2, and output files
    Read1 = input_folder + "/" + sample + ".R1.fastq.gz"
    Read2 = input_folder + "/" + sample + ".R3.fastq.gz"
    output_file = output_folder + "/" + sample + ".R2.fastq.gz"
    mismatch_rate = int(mismatch_rate)
    f1 = gzip.open(Read1, "rt")
    f2 = gzip.open(Read2, "rt")
    f3 = gzip.open(output_file, 'wt')

    line1 = f1.readline()
    line2 = f2.readline()

    total_line = 0
    filtered_line = 0

    while line1:
        total_line += 1
        line1 = f1.readline()

        bead_barcode_1 = line1[0:11]

        if bead_barcode_1 in list_barcode_1:
            bc_match_1 = list_barcode_1[bead_barcode_1]
            bc1_length = len(bc_match_1)

            bead_barcode_2 = line1[bc1_length + 4: bc1_length + 12]
            bead_barcode_3 = line1[bc1_length + 16: bc1_length + 24]
            bead_barcode_4 = line1[bc1_length + 34: bc1_length + 42]
            bead_UMI = line1[bc1_length + 42:]

            if (bead_barcode_2 in list_barcode_2) and (bead_barcode_3 in list_barcode_3) and (bead_barcode_4 in list_barcode_4):
                filtered_line += 1
                first_line = '@' + bc_match_1 + ',' + list_barcode_2[bead_barcode_2] + ',' + list_barcode_3[bead_barcode_3] + ','+ list_barcode_4[bead_barcode_4] + "," + bead_UMI + ',' + line2[1:]
                f3.write(first_line)
                
                second_line = f2.readline()
                f3.write(second_line)

                third_line = f2.readline()
                f3.write(third_line)

                four_line = f2.readline()
                f3.write(four_line)
                
                line2 = f2.readline()
                
            else:
                line2 = f2.readline()
                line2 = f2.readline()
                line2 = f2.readline()
                line2 = f2.readline()
        else:
            line2 = f2.readline()
            line2 = f2.readline()
            line2 = f2.readline()
            line2 = f2.readline()

        line1 = f1.readline()
        line1 = f1.readline()
        line1 = f1.readline()
        

        

    f1.close()
    f2.close()
    f3.close()
    print("sample name: %s, total line: %f, filtered line: %f, filter rate: %f"
          % (sample, total_line, filtered_line, float(filtered_line) / float(total_line)))


def extract_spatial_barcode_files(input_folder, sampleID, output_folder, core, list_barcode_1_file, list_barcode_2_file, list_barcode_3_file, list_barcode_4_file, mismatch_rate=1):
    init_message = '''
    --------------------------start attaching UMI-----------------------------
    input folder: %s
    sample ID: %s
    output_folder: %s
    Barcode 1 file: %s
    Barcode 2 file: %s
    Barcode 3 file: %s
    Barcode 4 file: %s
    ___________________________________________________________________________
    ''' % (input_folder, sampleID, output_folder, list_barcode_1_file, list_barcode_2_file, list_barcode_3_file, list_barcode_4_file)

    print(init_message)

    print("Load barcode dictionary...")
    # load the barcode dictionary
    with open(list_barcode_1_file, 'rb') as f:
        list_barcode_1 = pickle.load(f)

    with open(list_barcode_2_file, 'rb') as f:
        list_barcode_2 = pickle.load(f)

    with open(list_barcode_3_file, 'rb') as f:
        list_barcode_3 = pickle.load(f)

    with open(list_barcode_4_file, 'rb') as f:
        list_barcode_4 = pickle.load(f)

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
    func = partial(extract_spatial_barcode, input_folder = input_folder, output_folder=output_folder,
                   list_barcode_1=list_barcode_1, list_barcode_2=list_barcode_2,
                   list_barcode_3=list_barcode_3, list_barcode_4=list_barcode_4, mismatch_rate=mismatch_rate)
    result = p.map(func, sample_list)
    p.close()
    p.join()

    # print the completion message
    com_message = '''~~~~~~~~~~~~~~~Spatial barcode extraction done~~~~~~~~~~~~~~~~~~'''
    print(com_message)


if __name__ == "__main__":
    input_folder = sys.argv[1]
    sampleID = sys.argv[2]
    output_folder = sys.argv[3]
    core = sys.argv[4]
    list_barcode_1_file = sys.argv[5]
    list_barcode_2_file = sys.argv[6]
    list_barcode_3_file = sys.argv[7]
    list_barcode_4_file = sys.argv[8]
    extract_spatial_barcode_files(input_folder, sampleID, output_folder, core, list_barcode_1_file, list_barcode_2_file, list_barcode_3_file, list_barcode_4_file, mismatch_rate=1)

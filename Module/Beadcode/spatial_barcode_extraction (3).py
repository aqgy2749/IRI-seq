import subprocess
import sys
from Levenshtein import distance
import gzip
from multiprocessing import Pool
from functools import partial
import pickle

def extract_spatial_barcode(sample, input_folder, output_folder, list_barcode_1, list_barcode_2, list_barcode_3, list_barcode_4):
    # open the read1, read2, and output file
    Read1 = input_folder + "/" + sample + ".R2.fastq.gz"
    output_file = output_folder + "/" + sample + ".spatial.txt.gz"
    f1 = gzip.open(Read1, "rt")
    f3 = gzip.open(output_file, 'wt')

    total_line = 0
    filtered_line = 0

    while True:
        line1 = f1.readline()
        if not line1:  # If there are no more lines, exit the loop
            break

        total_line += 1

        # Process the first line
        Ligation = line1.strip().split(",")

        bc1 = Ligation[0][1:]
        bc2 = Ligation[1]
        bc3 = Ligation[2]
        bc4 = Ligation[3]
        UMI = Ligation[4]
        

        # Skip two lines
        f1.readline()

        # Read the third line
        line3 = f1.readline().strip()
        

        bead_barcode_1 = line3[0:11]

        if bead_barcode_1 in list_barcode_1:
            bc_match_1 = list_barcode_1[bead_barcode_1]
            bc1_length = len(bc_match_1)

            bead_barcode_2 = line3[bc1_length + 4: bc1_length + 12]
            bead_barcode_3 = line3[bc1_length + 16: bc1_length + 24]
            bead_barcode_4 = line3[bc1_length + 34: bc1_length + 46]
            

            if (bead_barcode_2 in list_barcode_2) and (bead_barcode_3 in list_barcode_3) and (bead_barcode_4 in list_barcode_4) :
                full_barcode = bc1 + bc2 +  bc3 + bc4 + "," + UMI + "," + bc_match_1  + list_barcode_2[bead_barcode_2] + list_barcode_3[bead_barcode_3]  + list_barcode_4[bead_barcode_4] 
                filtered_line += 1
                f3.write(full_barcode + "\n")
                
        f1.readline()
        f1.readline()

    f1.close()
    f3.close()

    print("sample name: %s, total line: %f, filtered line: %f, filter rate: %f"
          % (sample, total_line, filtered_line, float(filtered_line) / float(total_line)))

def extract_spatial_barcode_files(input_folder, sampleID, output_folder, core, list_barcode_1_file, list_barcode_2_file, list_barcode_3_file, list_barcode_4_file):
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

    with open(list_barcode_4_file2, 'rb') as f:
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
    func = partial(extract_spatial_barcode, input_folder=input_folder, output_folder=output_folder,
                   list_barcode_1 = list_barcode_1, list_barcode_2 = list_barcode_2, 
                   list_barcode_3 = list_barcode_3, list_barcode_4 = list_barcode_4)
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
    list_barcode_4_file2 = sys.argv[8]
    extract_spatial_barcode_files(input_folder, sampleID, output_folder, core, list_barcode_1_file, list_barcode_2_file, list_barcode_3_file, list_barcode_4_file2)

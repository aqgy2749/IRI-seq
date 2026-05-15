import os
import shutil
def dir_make(output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    else:
        print("Folder already exist!")
        
def dir_rm(output_folder):
    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)

def dir_cp(input_folder, output_folder):
    shutil.copytree(input_folder, output_folder)

def dir_mv(input_folder, output_folder):
    shutil.move(input_folder, output_folder)

def file_cp(input_file, output_file):
    shutil.copy(input_file, output_file)

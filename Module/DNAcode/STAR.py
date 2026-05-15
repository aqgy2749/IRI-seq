import subprocess
def Fastq_star_alignment_multi_files(input_folder, sample_ID, output_folder, core_num = 2,  
                         index = "/ru-auth/local/home/asziraki/projects/AS_20200820_single_cell_pipeline/Raw_data/AS_20200820_genome_files/original_pipeline/Reference_ALL/index/index/STAR/STAR_hg19_mm10_RNAseq/", 
                         star_path="/rugpfs/fs0/cao_lab/scratch/asziraki/anaconda3/envs/original_pipeline/bin/STAR"):
    input_command = f"for sample in $(cat {sample_ID}); do echo Aligning $sample; {star_path} --runThreadN {core_num} --outSAMstrandField intronMotif --genomeDir {index} --readFilesCommand zcat --readFilesIn {input_folder}/$sample*gz --outFileNamePrefix {output_folder}/$sample --genomeLoad LoadAndKeep; done"
    print(input_command)
    result = subprocess.check_output(input_command, shell=True, text=True)
    print(result)
    
    input_command = f"{star_path} --genomeDir {index} --genomeLoad Remove"
    print(input_command)
    result = subprocess.check_output(input_command, shell=True, text=True)
    print(result)

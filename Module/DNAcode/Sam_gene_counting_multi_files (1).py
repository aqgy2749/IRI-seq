import itertools
import collections
import os
import numpy as np
import pandas as pd
from multiprocessing import *
import gzip
import HTSeq
import sys
from functools import partial
import logging
import pickle

def find_intersected_gene(alnmt, exons, genes, gene_end, gene_start, exon_only):
    # Define a function that return the intersected gene for alignment read
    # each message include: chrom, start, end, read name, 0/1, strand, Exon/Intron information, Single/Multple genes, gene id, the id of gene with nearest 5' end, the id of gene with nearest 3' end, intersection exon ID, correct or opposite strand
    line_name = alnmt.get_sam_line().split("\t")[0]
    line_interval = alnmt.iv
    message = "\t".join(map(str, [line_interval.chrom, line_interval.start, line_interval.end, line_name, 1, line_interval.strand]))
    align_type = "no_match"
    
    if not alnmt.aligned:
        align_type = "no_match_aligned"
        message = "\t".join(map(str, ["Not_aligned", 0, 0, line_name, 0, "Not_aligned", "Not_aligned", "Not_aligned", "Not_aligned", "Not_aligned", "Not_aligned", "Not_aligned"]))
    elif alnmt.iv.chrom not in genes.chrom_vectors:
        align_type = "no_match_chrom"
        message = "\t".join(map(str,[line_interval.chrom, line_interval.start, line_interval.end, line_name, 1, line_interval.strand]))
        message = "\t".join(map(str,[message, "Not_aligned", "Not_aligned", "Not_aligned", "Not_aligned", "Not_aligned", "Not_aligned"]))
    else:
        # First check the intersectin with exons
        gene_id_intersect = set()
        gene_id_combine = set()
        inter_count = 0
        for cigop in alnmt.cigar:
            if cigop.type != "M":
                continue

            for iv,val in exons[cigop.ref_iv].steps():
                #print iv, val
                gene_id_combine |= val
                if inter_count == 0:
                    gene_id_intersect |= val
                    inter_count += 1
                else:
                    gene_id_intersect &= val
                #print "intersect set:", gene_id_intersect
                #print "combine set:", gene_id_combine
        # first check the intersection set
        if len(gene_id_intersect) == 1:
            gene_id = list(gene_id_intersect)[0]
            
            # find the aligned exon ID
            exon_id_intersect = set()
            exon_id_combine = set()
            inter_count_exon = 0
            
            for cigop in alnmt.cigar:
                if cigop.type != "M":
                    continue

                for iv,val in exon_only[cigop.ref_iv].steps():
                    #print iv, val
                    exon_id_combine |= val
                    if inter_count_exon == 0:
                        exon_id_intersect |= val
                        inter_count_exon += 1
                    else:
                        exon_id_intersect &= val
            
            if len(exon_id_intersect) == 1:
                exon_id = list(exon_id_intersect)[0]
                message = "\t".join(map(str,[message, "Exon", "Single_gene", gene_id, gene_id, gene_id, exon_id])) + "\n"
            else:
                message = "\t".join(map(str,[message, "Exon", "Single_gene", gene_id, gene_id, gene_id, "No_unique_exon_match"])) + "\n"
            align_type = "perfect_inter_exon"
            
        elif len(gene_id_intersect) > 1:
            gene_id_start = find_nearest_gene(alnmt.iv.start_d, gene_id_intersect, gene_start)
            gene_id_end = find_nearest_gene(alnmt.iv.end_d, gene_id_intersect, gene_end)
            message = "\t".join(map(str,[message, "Exon", "Multi_gene", ",".join(gene_id_intersect), gene_id_start, gene_id_end, "No_unique_exon_match"])) + "\n"
            align_type = "nearest_inter_exon"
        else:
            # if there no intersection match, then find the union sets
            if len(gene_id_combine) == 1:
                gene_id = list(gene_id_combine)[0]
                message = "\t".join(map(str,[message, "Exon", "Single_gene", gene_id, gene_id, gene_id, "No_unique_exon_match"])) + "\n"
                align_type = "perfect_combine_exon"
                
            elif len(gene_id_combine) > 1:
                gene_id_start = find_nearest_gene(alnmt.iv.start_d, gene_id_combine, gene_start)
                gene_id_end = find_nearest_gene(alnmt.iv.end_d, gene_id_combine, gene_end)
                message = "\t".join(map(str,[message, "Exon", "Multi_gene", ",".join(gene_id_combine), gene_id_start, gene_id_end, "No_unique_exon_match"])) + "\n"
                align_type = "nearest_combine_exon"
                
            else:
                # if there is no intersection match or union match, then search for genes to find the intronic match
                gene_id_intersect = set()
                gene_id_combine = set()
                inter_count = 0
                for cigop in alnmt.cigar:
                    if cigop.type != "M":
                        continue
                    for iv,val in genes[cigop.ref_iv].steps():
                        gene_id_combine |= val
                        if inter_count == 0:
                            gene_id_intersect |= val
                            inter_count += 1
                        else:
                            gene_id_intersect &= val

                if len(gene_id_intersect) == 1:
                    gene_id = list(gene_id_intersect)[0]
                    message = "\t".join([message, "Intron", "Single_gene", gene_id, gene_id, gene_id, "No_unique_exon_match"]) + "\n"
                    align_type = "perfect_inter_gene"

                elif len(gene_id_intersect) > 1:
                    gene_id_start = find_nearest_gene(alnmt.iv.start_d, gene_id_intersect, gene_start)
                    gene_id_end = find_nearest_gene(alnmt.iv.end_d, gene_id_intersect, gene_end)
                    message = "\t".join([message, "Intron", "Multi_gene", ",".join(gene_id_intersect), gene_id_start, gene_id_end, "No_unique_exon_match"]) + "\n"
                    align_type = "nearest_inter_gene"

                else:
                    # if there no intersection match, then find the union sets
                    if len(gene_id_combine) == 1:
                        gene_id = list(gene_id_combine)[0]
                        message = "\t".join([message, "Intron", "Single_gene", gene_id, gene_id, gene_id, "No_unique_exon_match"]) + "\n"
                        align_type = "perfect_combine_gene"

                    elif len(gene_id_combine) > 1:
                        gene_id_start = find_nearest_gene(alnmt.iv.start_d, gene_id_combine, gene_start)
                        gene_id_end = find_nearest_gene(alnmt.iv.end_d, gene_id_combine, gene_end)
                        message = "\t".join([message, "Intron", "Multi_gene", ",".join(gene_id_combine), gene_id_start, gene_id_end, "No_unique_exon_match"]) + "\n"
                        align_type = "nearest_combine_gene"

                    else:
                        align_type = "no_match_gene"
                        message = "\t".join(map(str,[message, "Not_aligned", "Not_aligned", "Not_aligned", "Not_aligned", "Not_aligned", "Not_aligned"]))
    return((message, align_type))



def scRNAseq_count(sample, input_folder, output_folder, exons, genes, gene_end, gene_start, exon_only):
    sam_file = input_folder + "/" + sample + ".sam"
    report = output_folder + "/" + sample + ".report"
    count_output = output_folder + "/" + sample + ".count.gz"
    counts = collections.Counter()
    total_count = 0
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    output_count = gzip.open(count_output, 'wt')
    
    print("Start read the input file: " + sam_file + "....")
    
    with HTSeq.SAM_Reader(sam_file) as almnt_file:
        for alnmt in almnt_file:
            total_count += 1
            aln_result = find_intersected_gene(alnmt, exons, genes, gene_end, gene_start, exon_only)
            if aln_result[1] != "no_match_gene":
                counts[aln_result[1]] += 1
                output_count.write(aln_result[0].strip() + "\tCorrect_strand\n")
            else:
                if alnmt.iv.strand == "+":
                    alnmt.iv.strand = "-"
                else:
                    alnmt.iv.strand = "+"
                
                aln_result_oppo_strand = find_intersected_gene(alnmt, exons, genes, gene_end, gene_start, exon_only)
                
                if aln_result_oppo_strand[1] == "no_match_gene":
                    counts[aln_result[1]] += 1
                    output_count.write(aln_result[0] + "\tNo_match\n")
                else:
                    counts[aln_result_oppo_strand[1]] += 1
                    
                    message = "\t".join([aln_result_oppo_strand[0].strip(), "Opposite_strand\n"])
                    output_count.write(message)
                    
    output_count.close()
    
    print("File name: ", sam_file)
    print("1: Perfect intersect exon match: ", counts["perfect_inter_exon"])
    print("2: Nearest intersect exon match: ", counts["nearest_inter_exon"])
    print("3: Perfect combine exon match: ", counts["perfect_combine_exon"])
    print("4: Nearest combine exon match: ", counts["nearest_combine_exon"])
    print("5: Perfect intersect gene match: ", counts["perfect_inter_gene"])
    print("6: Nearest intersect gene match: ", counts["nearest_inter_gene"])
    print("7: Perfect combine gene match: ", counts["perfect_combine_gene"])
    print("8: Nearest combine gene match: ", counts["nearest_combine_gene"])
    print("9: No alignment: ",  counts["no_match_aligned"])
    print("10: Wrong chromosome: ",  counts["no_match_chrom"])
    print("11: No match gene: ",  counts["no_match_gene"])
    print("12: Total count: ", total_count)
    print("Sam file analysis finished~")
    
    with open(report, 'w') as report:
        report.write("1: Perfect intersect exon match: " + "," + str(sample) + "," + str(counts["perfect_inter_exon"]) + "\n")
        report.write("2: Nearest intersect exon match: " + "," + str(sample) + "," + str(counts["nearest_inter_exon"]) + "\n")
        report.write("3: Perfect combine exon match: " + "," + str(sample) + "," + str(counts["perfect_combine_exon"]) + "\n")
        report.write("4: Nearest combine exon match: " + "," + str(sample) + "," + str(counts["nearest_combine_exon"]) + "\n")
        report.write("5: Perfect intersect gene match: " + "," + str(sample) + "," + str(counts["perfect_inter_gene"]) + "\n")
        report.write("6: Nearest intersect gene match: " + "," + str(sample) + "," + str(counts["nearest_inter_gene"]) + "\n")
        report.write("7: Perfect combine gene match: " + "," + str(sample) + "," + str(counts["perfect_combine_gene"]) + "\n")
        report.write("8: Nearest combine gene match: " + "," + str(sample) + "," + str(counts["nearest_combine_gene"]) + "\n")
        report.write("9: No match: " + "," + str(sample) + "," + str(counts["no_match_aligned"]) + "\n")
        report.write("10: Wrong chromosome: " + "," + str(sample) + "," + str(counts["no_match_chrom"]) + "\n")
        report.write("11: No match gene: " + "," + str(sample) + "," + str(counts["no_match_gene"]) + "\n")
        report.write("12: Total count: " + "," + str(sample) + "," + str(total_count) + "\n")
    
    return 0

def find_nearest_gene(al_end, gene_id_intersect, gene_end):
    gene_id_end = {}
    for gene in gene_id_intersect:
        if gene in gene_end:
            gene_id_end[gene] = (abs(np.array(list(gene_end[gene])) - al_end)).min()
        else:
            print("****************Found one gene without transcript annotation*****************", "Gene name: ", gene)
    # filter the gene with the least distance. If there are two genes with the least distance, then  "_ambiguous" 
    # would be returned
    min_distance = min(gene_id_end.values())
    gene_ids_with_min_distance = []
    for gene_id, distance in gene_id_end.items():
        if distance == min_distance:
            gene_ids_with_min_distance.append(gene_id)
    
    if len(gene_ids_with_min_distance) > 1:
        gene_id = "_ambiguous"
    else:
        gene_id = gene_ids_with_min_distance[0]
    
    return(gene_id)

def scRNA_count_parallel(input_folder, sample_ID, output_folder, reference_file, core_number):
    # read in the gtf file, and then construct the genome interval for exons, genes, and gene end dictionary
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    print("Load gene reference file....")
    with open(reference_file, "rb") as f1:
        gene_reference = pickle.load(f1)
    genes = gene_reference["genes"]
    exons = gene_reference["exons"]
    exon_only = gene_reference["exon_only"]
    gene_end = gene_reference["gene_end"]
    gene_start = gene_reference["gene_start"]
    gene_annotat =  gene_reference["gene_annotat"]
    exon_annotate = gene_reference["exon_annotat"]
    
    print("Start processing all file....")
    sample_ID = list(pd.read_csv(sample_ID, header=None)[0])
    
    #scRNAseq_count(sample_ID[0], input_folder = input_folder, output_folder = output_folder, exons=exons, genes=genes, gene_end = gene_end, gene_start = gene_start, exon_only = exon_only)

    # parallele for the functions
    p = Pool(processes = int(core_number))
    #print("Processing core number: ", core_number)
    func = partial(scRNAseq_count, input_folder = input_folder, output_folder = output_folder, exons=exons, genes=genes, gene_end = gene_end, gene_start = gene_start, exon_only = exon_only)
    #sciRNAseq_count(sample, input_folder, exons, genes, gene_end, exon_only)
    result = p.map(func, sample_ID)
    p.close()
    p.join()
    print("Write gene annotation file into the output folder....")
    gene_annotat.to_csv(os.path.join(output_folder, "gene_anno.csv"), index = False)
    
    print("Write exon annotation file into the output folder....")
    exon_annotate.to_csv(os.path.join(output_folder, "exon_anno.csv"), index = False)
    
    print("All analysis done~")
'''    
if __name__ == "__main__":
    input_folder = sys.argv[1]
    sample_ID = sys.argv[2]
    output_folder = sys.argv[3]
    reference_file = sys.argv[4]
    core_number = sys.argv[5]
    scRNA_count_parallel(input_folder, sample_ID, output_folder, reference_file, core_number)
'''
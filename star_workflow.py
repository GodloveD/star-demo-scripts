#!/usr/bin/env python3

import subprocess
import os
import re
import glob
import sys
import datetime
from multiprocessing import Pool
from collections import defaultdict



# def create_results_dirs(fastqs_dir, results_dir):
#     os.makedirs(results_dir, exist_ok = True)
#     sample_re = re.compile("([^/]+).fastq.gz$")
#     input_file_list = set([sample_re.search(file).group(1) for file in glob.glob(os.path.join(fastqs_dir, "*.fastq.gz"))])
# 
#     # debugging
#     print(input_file_list)
# 
#     for file in input_file_list:
#         results_dir_path = os.path.join(results_dir, file)
#         os.makedirs(results_dir_path, exist_ok = True)

# loads the genome into shared memory so all STAR procs can use that copy
def startup(genome_dir, wdir):
    # mode here is important - "LoadAndExit" specifies that it should be stored 
    # in shared mem, would be "NoSharedMemory" otherwise
    print("Loading genome...")
    print('subprocess.run(["STAR", "--genomeDir", ' + genome_dir + ', "--genomeLoad", "LoadAndExit"], cwd=' + wdir + ', shell=False)')
    subprocess.run(["STAR", "--genomeDir", genome_dir, "--genomeLoad", "LoadAndExit"], cwd=wdir, shell=False)
    print("Genome loaded!")
    
# takes a sample folder name
def run_workflow(file):

    velocyto = False
    h5ad = True

    print("Starting workflow...")

    # get all the fastqs from the fastq storage corresponding to the folder
    # name given
    # sample_files = [f for f in glob.glob("/data/fastqs/" + results_dir + "*")]
    print(file)
    sample_file = os.path.join(fastqs_dir, fname  + ".fastq.gz")
    print(sample_file)

    # make a subdirectory for results
    # this is where we should make the first results dir
    cwd_path = os.path.join(results_dir, fname, "Pass1")
    os.makedirs(cwd_path, exist_ok = True)
 
    print("Running STAR...")

    # runs STAR - note that "LoadAndKeep" signals that the shared memory genome 
    # should be used, and that you still have to specify the same genome 
    # location
    # the number of procs used for the alignment of the sample folder's files 
    # can be set at "--runThreadN"
    # "readFilesIn" needs a comma with no space between the two fastq file 
    # paths
    # "--readFilesCommand" is currently set to indicate fastqs will be inside 
    # gzipped folders - STAR will do decompression so we don't need to store 
    # them uncompressed
    cmd = ["STAR",
           "--outFilterType",
           "BySJout",
           "--outFilterMultimapNmax",
           "20",
           "--alignSJoverhangMin",
           "8",
           "--alignSJDBoverhangMin",
           "1",
           "--outFilterMismatchNmax",
           "999",
           "--outFilterMismatchNoverLmax",
           "0.04",
           "--alignIntronMin",
           "20",
           "--alignIntronMax",
           "1000000",
           "--alignMatesGapMax",
           "1000000",
           "--outSAMstrandField",
           "intronMotif",
           "--outSAMtype",
           "BAM",
           "Unsorted",
           "--outSAMattributes",
           "NH",
           "HI",
           "NM",
           "MD",
           "--genomeLoad",
           "LoadAndKeep",
           "--outReadsUnmapped",
           "Fastx",
           "--readFilesCommand",
           "zcat",
           "--runThreadN",
           str(procs),
           "--genomeDir",
           genome_dir,
           "--readFilesIn",
           sample_file
           ]

    print(cmd)
    subprocess.run(cmd, cwd=cwd_path, shell=False)
# 
#     print("Running samtools sort...")
#     
#     #intermediate samtools steps that run sequentially and only take seconds to complete for each sample
#     subprocess.run(["samtools", 
#                     "sort", 
#                     "-m", 
#                     "6000000000", 
#                     "-o", 
#                     "./Aligned.out.sorted.bam",
#                     "./Aligned.out.bam"
#                     ], cwd=cwd_path, shell=False)
# 
#     print("Running samtools index...")
 
#     subprocess.run(["samtools", 
#                     "index",
#                     "-b",
#                     "./Aligned.out.sorted.bam"
#                     ], cwd=cwd_path, shell=False)
# 
# 
#     if(h5ad):
#         print("Running samtools htgen...")
#     
#         subprocess.run(["samtools",
#                         "sort",
#                         "-m",
#                         "6000000000",
#                         "-n",
#                         "-o",
#                         "./Aligned.out.sorted-byname.bam",
#                         "./Aligned.out.sorted.bam"
#                         ], cwd=cwd_path, shell=False)
# 
#         htseq_out = open("/data/results/" + results_dir + "/Pass1/htseq-count.txt", "w+")
# 
#         print("Running htseq...")
# 
#         #htseq generates a file "htseq-count.txt" with the genes defined in the GTF (I think) and the number of times they were found expressed in the sample
#         #this is the final product and what the h5ad file is made from
#         subprocess.run(["htseq-count", 
#                         "-r",
#                         "name",
#                         "-s"
#                         "no",
#                         "-f",
#                         "bam",
#                         "--idattr=gene_id",
#                         "-m",
#                         "intersection-nonempty",
#                         "./Aligned.out.sorted-byname.bam",
#                         "/data/STAR/Saccharomyces_cerevisiae.R64-1-1.79.gtf",
#                         ], cwd=cwd_path, shell=False, stdout=htseq_out)
#     if(velocyto):
#         bam_files = [f for f in glob.glob("/data/results/" + results_dir + "/Pass1/Aligned.out.sorted.bam")]
#         #add velocyto run here - only works over human and mouse genomes
#         subprocess.run()
# 
#     print("Done!")


if __name__ == "__main__":                                                                                            
    global fastqs_dir 
    global results_dir
    global genome_dir
    global procs
    fastqs_dir = "/data/fastqs"
    genome_dir = "/data/STAR/genome"
    results_dir = "/data/results"                                                      
    wdir = "/data"
    procs = 8
    # if len(sys.argv) > 5:                                                       
    #     print("ERROR: too many args")                                           
    # elif len(sys.argv) == 5:
    #     procs = sys.argv[4]
    # elif len(sys.argv) == 4:
    #     wdir = syst.argv[3]
    # elif len(sys.argv) == 3:                                                    
    #     results_dir = sys.argv[2]                                                      
    # elif len(sys.argv) == 2:                                                    
    #     fastqs_dir = sys.argv[1]                                                      
                                                                                
    # create_results_dirs(fastqs_dir, results_dir) 

    #load genome
    startup(genome_dir, wdir)

    # create processor pool
    pool = Pool(processes=int(procs))
    print("Processor pool loaded!")

    #grab all the input folders that were created by the preprocessing step
    # input_folder_list = sorted([os.path.basename(f) for f in glob.glob(os.path.join(results_dir, "*"))])
    # input_folder_list = sorted(glob.glob(os.path.join(results_dir, "*")))
    sample_re = re.compile("([^/]+).fastq.gz$")
    input_file_list = set([sample_re.search(file).group(1) for file in glob.glob(os.path.join(fastqs_dir, "*.fastq.gz"))])
    print(input_file_list)
    print("Inputs enumerated with len: " + str(len(input_file_list)))
    pool.map(run_workflow, input_file_list)
    #run the workflow - pool.map basically says run this method (run_workflow) with this set of inputs (input_folder_list) using the processors we defined (16)
    #might need a couple vars here passed in from the workflow to tell whether the h5ad, loom, or both types of output files should be generated

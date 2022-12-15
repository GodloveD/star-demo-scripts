#!/usr/bin/env python3

import subprocess
import os
import re
import glob
import argparse
from multiprocessing import Pool
from collections import defaultdict


# loads the genome into shared memory so all STAR procs can use that copy
def startup():
    # mode here is important - "LoadAndExit" specifies that it should be stored 
    # in shared mem, would be "NoSharedMemory" otherwise
    print("Loading genome...")
    cmd = ["STAR", "--genomeDir", genome_dir, "--genomeLoad", "LoadAndExit"]
    # print(cmd)
    subprocess.run(cmd, cwd=wdir, shell=False)
    print("Genome loaded!")
    

# takes a sample folder name
def run_workflow(file):

    print("Starting workflow...")

    # get all the fastqs from the fastq storage corresponding to the folder
    # name given
    # sample_files = [f for f in glob.glob("/data/fastqs/" + results_dir + "*")]
    # print(file)
    sample_file = os.path.join(fastqs_dir, file  + ".fastq.gz")
    # print(sample_file)

    # make a subdirectory for results
    # this is where we should make the first results dir
    cwd_path = os.path.join(results_dir, file, "Pass1")
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
           str(threads),
           "--genomeDir",
           genome_dir,
           "--readFilesIn",
           sample_file
           ]

    # print(cmd)
    subprocess.run(cmd, cwd=cwd_path, shell=False)

    print("Running samtools sort...")
    
    # intermediate samtools steps that run sequentially and only take seconds 
    # to complete for each sample
    subprocess.run(["samtools", 
                    "sort", 
                    "-m", 
                    "6000000000", 
                    "-o", 
                    "./Aligned.out.sorted.bam",
                    "./Aligned.out.bam"
                    ], cwd=cwd_path, shell=False)

    print("Running samtools index...")

    subprocess.run(["samtools", 
                    "index",
                    "-b",
                    "./Aligned.out.sorted.bam"
                    ], cwd=cwd_path, shell=False)

    print("Running samtools htgen...")
    
    subprocess.run(["samtools",
                    "sort",
                    "-m",
                    "6000000000",
                    "-n",
                    "-o",
                    "./Aligned.out.sorted-byname.bam",
                    "./Aligned.out.sorted.bam"
                    ], cwd=cwd_path, shell=False)

    htseq_out = open(os.path.join(cwd_path, "htseq-count.txt"), "w+")

    print("Running htseq...")

    # htseq generates a file "htseq-count.txt" with the genes defined in the GTF 
    # (I think) and the number of times they were found expressed in the sample
    # this is the final product and what the h5ad file is made from
    subprocess.run(["htseq-count", 
                        "-r",
                        "name",
                        "-s"
                        "no",
                        "-f",
                        "bam",
                        "--idattr=gene_id",
                        "-m",
                        "intersection-nonempty",
                        "./Aligned.out.sorted-byname.bam",
                        os.path.join(ref_dir, "Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.21.gtf"),
                        ], cwd=cwd_path, shell=False, stdout=htseq_out)


if __name__ == "__main__":  

    global wdir
    global fastqs_dir 
    global ref_dir
    global genome_dir
    global results_dir
    global threads
    global gtf_file
    
    parser = argparse.ArgumentParser()

    parser.add_argument('--working_dir',   
        type=str, 
        help='''the parent directory where the data are located (in 
             subdirectories) and where the analysis will be carried out''', 
        default="/data")
    parser.add_argument('--fastqs_dir',    
        type=str, 
        help='''the directory where the fastq samples have been downloaded''', 
        default="/data/fastqs")
    parser.add_argument('--reference_dir', 
        type=str, 
        help='''the directory where the raw reference genome has be downloaded 
             and gunzipped''', 
        default="/data/ref")
    parser.add_argument('--genome_dir',    
        type=str, 
        help='''the location of the STAR indexed genome''', 
        default="/data/STAR/genome")
    parser.add_argument('--results_dir',   
        type=str, 
        help='''where to put the data analyzed by this script''', 
        default="/data/results")
    parser.add_argument('--procs',         
        type=int, 
        help='''the number of processes (for simultaneous STAR runs sharing the 
             loaded reference) defaults to the number of files''')
    parser.add_argument('--threads',         
        type=int, 
        help='''the number of threads for each STAR process''', 
        default=4)
    parser.add_argument('--gtf_file',         
        type=str, 
        help='''the name of the gtf index file (assumed to be in the ref_dir))''', 
        default="Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.21.gtf")

    args = parser.parse_args()

    wdir        = args.working_dir
    fastqs_dir  = args.fastqs_dir
    ref_dir     = args.reference_dir
    genome_dir  = args.genome_dir
    results_dir = args.results_dir
    procs       = args.procs
    threads     = args.threads
    gtf_file    = args.gtf_file

    #load genome
    startup()

    # grab all the input files that were downloaded
    sample_re = re.compile("([^/]+).fastq.gz$")
    input_file_list = set([sample_re.search(file).group(1) for file in glob.glob(os.path.join(fastqs_dir, "*.fastq.gz"))])
    # print(input_file_list)
    print("Inputs enumerated with len: " + str(len(input_file_list)))

    # create processor pool
    if procs is None:
        pool = Pool(processes=len(input_file_list))
    else:
        pool = Pool(processes=procs)
    print("Processor pool loaded!")
    
    # run the workflow - pool.map basically says run this method (run_workflow) 
    # with this set of inputs (input_folder_list) using the processors we 
    # defined (default 16)
    pool.map(run_workflow, input_file_list)
    
    print("Done!")
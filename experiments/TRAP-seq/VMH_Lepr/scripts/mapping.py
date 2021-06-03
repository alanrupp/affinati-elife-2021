#!/usr/bin/python

'''
Analyzing FASTQ files to generate count tables. Running QC using FastQC,
filtering low quality reads with fastq_quality_filter, and mapping and counting
reads with STAR.
'''

import subprocess
import glob
import datetime
import re
import os
import pandas as pd
import numpy as np

# mkdir function
def mkdir(directory):
    out = f"mkdir {directory}"
    subprocess.call(out, shell=True)
    return out

# grab software version
def grab_version(program):
    arg = f"{program} --version"
    output = subprocess.Popen(arg, stdout=subprocess.PIPE, shell=True).communicate()[0]
    output = output.decode('utf-8')
    return output

# FASTQC
def fastqc(dir, sample, read):
    log = grab_version("fastqc")
    # pipe files through zcat and FASTQC into a subfolder called FASTQC / sample
    if not os.path.isdir(f"FASTQC/{sample}"):
        log += mkdir(f"FASTQC/{sample}") + "\n"
    log += mkdir(f"FASTQC/{sample}/{read}") + "\n"
    fastq = glob.glob(f"{dir}/{sample}/*{read}*fastq.gz")
    fastq = ' '.join(fastq)
    out = f"zcat {fastq} | "
    out += f"fastqc -t 8 -o FASTQC/{sample}/{read} stdin"
    # grab date and time
    log += f"{datetime.datetime.now().strftime('%c')}\n"
    print(f"{log}\n{out}\n")
    subprocess.call(out, shell=True)
    log += f"{out}\n"
    return f"{log}\n"

# FASTQ quality filter
def filter(dir, sample, read):
    log = 'fastq_quality_filter version 0.0.14\n'
    fastq = glob.glob(f"{dir}/{sample}/*{read}*fastq.gz")
    fastq = ' '.join(fastq)
    # pipe files through zcat and fastq_quality_filter
    out = f"zcat {fastq} | "
    out += "fastq_quality_filter -q 20 -z "
    out += f"-o {dir}/{sample}/{read}.fastq.gz"
    # grab date and time
    log += f"{datetime.datetime.now().strftime('%c')}\n"
    print(f"{log}\n{out}\n")
    subprocess.call(out, shell=True)
    log += f"{out}\n"
    return f"{log}\n"

# make a FASTA with custom genes
def custom_fasta(fasta_files):
    out = f"cat {' '.join(fasta_files)} > data/genome.fa"
    print(out)
    subprocess.call(out, shell=True)
    return f"{out}\n"

# make a GTF with custom genes
def custom_gtf(gtf_files):
    out = f"cat {' '.join(gtf_files)} > data/genes.gtf"
    print(out)
    subprocess.call(out, shell=True)
    return f"{out}\n"

# STAR genomeGenerate
def generate_genome(fasta, gtf):
    star = '/opt/STAR-2.7.3a/bin/Linux_x86_64_static/STAR'
    log = grab_version(star)
    if not os.path.isdir("data/genome"):
        log += f"{mkdir('data/genome')}\n"
    out = f"{star} "
    out += "--runMode genomeGenerate "
    out += "--genomeDir data/genome "
    out += "--runThreadN 8 "
    out += f"--genomeFastaFiles {fasta} "
    out += f"--sjdbGTFfile {gtf} "
    log += f"{datetime.datetime.now().strftime('%c')}\n"
    print(f"{log}\n{out}\n")
    subprocess.call(out, shell=True)
    log += f"{out}\n"
    return f"{log}\n"

# STAR GeneCounts
def count(dir, sample, gtf, paired_end=False):
    star = '/opt/STAR-2.7.3a/bin/Linux_x86_64_static/STAR'
    log = f"STAR {grab_version(star)}"
    out = f"{star} "
    out += "--runThreadN 8 "
    out += "--genomeDir data/genome "
    out += f"--sjdbGTFfile {gtf} "
    out += "--readFilesCommand zcat "
    if paired_end:
        out += f"--readFilesIn {dir}/{sample}/R1.fastq.gz {dir}/{sample}/R2.fastq.gz "
    else:
        out += f"--readFilesIn {dir}/{sample}/R1.fastq.gz "
    out += f"--outFileNamePrefix data/STAR/{sample} "
    out += "--outSAMtype None "
    out += "--quantMode GeneCounts"
    # grab date and time
    log += f"{datetime.datetime.now().strftime('%c')}\n"
    print(f"{log}\n{out}\n")
    subprocess.call(out, shell=True)
    log += f"{out}\n"
    return f"{log}\n"

def gunzip(file):
    out = f"gunzip {file}"
    print(out)
    subprocess.call(out, shell=True)
    return f"{out}\n"

def gzip(file):
    out = f"gzip {file}"
    print(out)
    subprocess.call(out, shell=True)
    return f"{out}\n"

def tx2gene(gtf):
    log = f"\n----- Generating tx2gene CSV file from {gtf} -----\n"
    fields = ['transcript_id', 'transcript_name', 'gene_id', 'gene_name']
    transcript_id = np.where(np.array(fields) == 'transcript_id')[0][0]
    # function to grab info associated with user input fields
    def grab_info(field):
        info = re.findall(field + ' \"(.*?)\"', line)
        if len(info) == 0:
            info = ""
        else:
            info = info[0]
        return info
    # generate dictionary from GTF file
    d = dict()
    if re.search("gz$", gtf):
        log += gunzip(gtf)
    with open(gtf, 'r') as f:
        for line in f:
            if line.startswith('#'): continue # skip header lines
            line = line.split(sep='\t')[8] # line with relevant info
            field_stash = list()
            for field in fields:
                hit = grab_info(field)
                field_stash.append(hit)
            # if the transcript_id is new, add it to the dictionary with its info
            if field_stash[transcript_id] not in d.keys():
                d[field_stash[transcript_id]] = np.array(field_stash)[np.arange(len(field_stash)) != transcript_id]
            else: continue
    if re.search("gz$", gtf):
        log += gzip(re.sub(".gz$", "", gtf))
    # turn into data frame for easy output
    df = pd.DataFrame.from_dict(d, orient='index')
    df.reset_index(level=0, inplace=True)
    rename_col = {'index': 'transcript_id'}
    count = 0
    for field in fields:
        if field == 'transcript_id': continue
        else:
            rename_col[count] = field
            count += 1
    df.rename(columns=rename_col, inplace=True)
    # write to CSV file
    log += "Saving file as data/tx2gene.csv"
    df.to_csv("data/tx2gene.csv", index=False)
    return f"{log}\n"

# move a file
def mv(file, destination):
    out = f"mv {file} {destination}"
    print(out)
    subprocess.call(out, shell=True)
    return f"{out}\n"

# remove a file
def rm(dir, recursive=False):
    flag = '-r' if recursive else ''
    out = f"rm {flag} {dir}"
    print(out)
    subprocess.call(out, shell=True)
    return f"{out}\n"

# - run -----------------------------------------------------------------------
if __name__ == '__main__':
    # set up command line arguments
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--directory', help='directory')
    parser.add_argument('--genome', help='genome')
    parser.add_argument('--gtf', help='gtf')
    parser.add_argument('--paired_end', action='store_true')
    parser.add_argument('--custom_reference', action='store_true')
    parser.add_argument('--clean', action='store_true')
    args = parser.parse_args()
    # read in run information and reformat
    dir = args.directory
    samples = glob.glob(f"{dir}/*")
    samples = [re.findall("Sample.+", x)[0] for x in samples]
    # initialize log file
    log = '-- Starting pipeline --\n'
    log += f"{datetime.datetime.now().strftime('%c')}\n"
    # make output folders
    log += f"{mkdir('FASTQC')}\n"
    # run FASTQC and fastq_quality_filter on all FASTQ files
    reads = ['R1', 'R2'] if args.paired_end else ['R1']
    for sample in samples:
        for read in reads:
            log += f"\n------ QC and filtering {sample} {read} ------\n"
            log += fastqc(dir, sample, read)
            log += filter(dir, sample, read)
    if args.custom_reference:
        # make custom FASTA and GTF with Cre & GFP:L10a
        log += "\n----- Making custom FASTA & GTF -----\n"
        fasta = glob.glob("data/*.fa*")
        # if gzipped, unzip them and change FASTA list to open with `cat`
        if any([re.search("gz$", x) for x in fasta]):
            gzipped = []
            for fa in fasta:
                if re.search("gz$", fa):
                    gzipped.append(fa)
                    fasta.remove(fa)
                    fasta.append(re.sub(".gz$", "", fa))
                    log += gunzip(fa)
            log += custom_fasta(fasta)
            # re-gzip the files that were originally gzipped
            for gz in gzipped:
                log += gzip(re.sub(".gz$", "", gz))
        # do the same for GTF
        gtf = glob.glob("data/*.gtf*")
        # if gzipped, unzip them and change GTF list to open with `cat`
        if any([re.search("gz$", x) for x in gtf]):
            gzipped = []
            for g in gtf:
                if re.search("gz$", g):
                    gzipped.append(g)
                    gtf.remove(g)
                    gtf.append(re.sub(".gz$", "", g))
                    log += gunzip(g)
            log += custom_gtf(gtf)
            # re gzip the files that were originally gzipped
            for gz in gzipped:
                log += gzip(re.sub(".gz$", "", gz))
        # generate genome with custom files
        log += "\n----- Generating genome index -----\n"
        log += generate_genome("data/genome.fa", "data/genes.gtf")
    else:
        # generate genome
        log += "\n----- Generating genome index -----\n"
        log += generate_genome(args.genome, args.gtf)
    # map reads to genome index and count by gene
    log += f"{mkdir('data/STAR')}\n"
    for sample in samples:
        log += f"\n------ Counting {sample} -----\n"
        if args.custom_reference:
            log += count(dir, sample, "data/genes.gtf", args.paired_end)
        else:
            log += count(dir, sample, args.gtf, args.paired_end)
    if args.clean:
        counts = glob.glob(f"data/STAR/*ReadsPerGene*")
        # move count files to new directory
        if len(counts) == len(samples):
            log += mkdir("data/counts")
            for count in counts:
                log += mv(count)
        # make tx2gene file for mapping gene_id and gene_name
        if args.custom_reference:
            log += tx2gene("data/genes.gtf")
        else:
            log += tx2gene(args.gtf)
        # delete other unnecessary files
        log += "\n----- Deleting unnecessary files -----\n"
        log += rm("data/STAR", recursive=True)
        log += rm("data/genome", recursive=True)
        #log += rm("_STARtmp", recursive=True)
        filtered_fastq = glob.glob("data/*/*/Sample*/R[1-2].fastq.gz")
        for f in filtered_fastq:
            log += rm(f)
        if args.custom_reference:
            log += rm("data/genome.fa")
            log += rm("data/genes.gtf")
    # write log file
    with open('scripts/log.txt', 'w') as f:
        f.write(log)

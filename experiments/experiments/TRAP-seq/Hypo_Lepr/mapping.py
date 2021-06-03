#!/usr/bin/python

"""
Analyzing data from Allison et al manuscript. This script maps the raw FASTQ
files to the mouse genome and generates count tables. All other analysis was
done in R.

FASTQ files are in nested folders with `FASTQ / Run / Sample / *fastq.gz` format

First run FastQC, then fastq_quality_filter, then STAR
"""

import subprocess
import glob
import datetime

# mkdir function
def mkdir(directory):
    out = 'mkdir ' + directory
    subprocess.call(out, shell=True)
    return(out)

# grab software version
def grab_version(program):
    arg = program + ' --version'
    output = subprocess.Popen(arg, stdout=subprocess.PIPE, shell=True).communicate()[0]
    output = output.decode('utf-8')
    return(output)

# FASTQC
def fastqc(run, sample):
    log = grab_version("fastqc")
    # pipe files through zcat and FASTQC into a subfolder called FASTQC / sample
    log += mkdir("FASTQC/" + sample) + '\n'
    out = "zcat " + 'FASTQ/' + run + '/' + sample + "/*.fastq.gz"
    out += " | fastqc -t 8 -o " + "FASTQC/" + sample + " stdin"
    # grab date and time
    log += datetime.datetime.now().strftime("%c") + '\n'
    print(log + '\n' + out + '\n')
    subprocess.call(out, shell=True)
    log += out + '\n'
    return(log + '\n')

# FASTQ quality filter
def filter(run, sample):
    log = 'fastq_quality_filter version 0.0.14\n'
    fastq = glob.glob('FASTQ/' + run + '/' + sample + '/' + '*fastq.gz')
    fastq = ' '.join(fastq)
    # pipe files through zcat and fastq_quality_filter
    out = "zcat " + fastq
    out += " | fastq_quality_filter -q 20 -z"
    out += " -o " + "FASTQ/" + run + "/" + sample + "/filtered.fastq.gz"
    # grab date and time
    log += datetime.datetime.now().strftime("%c") + '\n'
    print(log + '\n' + out + '\n')
    subprocess.call(out, shell=True)
    log += out + '\n'
    return(log + '\n')

# STAR alignment
def STAR(run, sample):
    star = '~/STAR/bin/Linux_x86_64/STAR'
    log = grab_version(star)
    log += mkdir('STAR_outs/' + sample) + '\n'
    out = star
    out += " --runThreadN 8"
    out += " --genomeDir /home/alanrupp/STAR/GRCm38-92_Cre_GfpL10a"
    out += " --sjdbGTFfile /home/alanrupp/STAR/GRCm38.92_Cre_GfpL10a.gtf"
    out += " --readFilesCommand zcat"
    out += " --readFilesIn " + 'FASTQ/' + run + '/' + sample + '/filtered.fastq.gz'
    out += " --outFileNamePrefix STAR_outs/" + sample
    out += " --outSAMtype BAM SortedByCoordinate"
    out += " --quantMode GeneCounts"
    out += " --twopassMode Basic"
    # grab date and time
    log += datetime.datetime.now().strftime("%c") + '\n'
    print(log + '\n' + out + '\n')
    subprocess.call(out, shell=True)
    log += out + '\n'
    return(log + '\n')


# - run -----------------------------------------------------------------------
if __name__ == '__main__':
    # read in run information
    import pandas as pd
    runs = pd.read_csv('runs.csv')

    # initialize log file
    log = 'Starting pipeline\n'
    log += datetime.datetime.now().strftime("%c") + '\n'

    # make output folders
    log += mkdir("FASTQC") + '\n'
    log += mkdir("STAR_outs") + '\n'

    # cycle through all samples
    for i, df in runs.iterrows():
        log += '\n------ ' + df['Sample'] + ' ------\n'
        log += fastqc(df['Run'], df['Sample'])
        log += filter(df['Run'], df['Sample'])
        log += STAR(df['Run'], df['Sample'])

    # write log file
    with open('mapping_log.txt', 'w') as f:
        f.write(log)

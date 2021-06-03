"""
Reorganizing the output files that result from mapping. I only want to keep the
count files from STAR and the original FASTQ files
"""

import subprocess
import glob
import datetime

# mkdir function
def mkdir(directory):
    out = 'mkdir ' + directory
    subprocess.call(out, shell=True)
    return(out + '\n')

# move files
def mv(file, directory):
    out = 'mv ' + file + ' ' + directory
    subprocess.call(out, shell=True)
    return(out + '\n')

# delete filtered fastq and other STAR outputs
def scrub(sample):
    star_files = glob.glob('STAR_outs/' + sample + '*')[0]
    star = 'rm -r ' + star_files
    #subprocess.call(star, shell=True)
    filtered_file = glob.glob('FASTQ/' + sample + '/filtered.fastq.gz')[0]
    filtered = 'rm ' + filtered_file
    #suprocess.call(filtered, shell=True)
    return(star + '\n' + filtered + '\n')

# - Main ----------------------------------------------------------------------
if __name__ == '__main__':
    import re
    log = 'Reorganizing files\n'
    log += datetime.datetime.now().strftime("%c") + '\n'

    # move count files into a new folder
    log += mkdir("Counts") + '\n'
    count_files = glob.glob('STAR_outs/*ReadsPerGene.out.tab')
    for file in count_files:
        log += mv(file, "Counts")

    # delete other associated files if counts files was successful
    counts = glob.glob("Counts/*")
    samples = [re.findall("(?<=Counts\\/)(.+)(?=ReadsPerGene)", x)[0] for x in counts]
    for sample in samples:
        log += scrub(sample)

    # write log file
    with open('organization_log.txt', 'w') as f:
        f.write(log)

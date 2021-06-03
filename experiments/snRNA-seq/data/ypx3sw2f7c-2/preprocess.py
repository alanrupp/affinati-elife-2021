#!/usr/bin/python

'''
Generating a metadata file from the data in this folder and gzipping the files
to save space
'''

import glob
import subprocess
import pandas as pd

def view_archive(file):
    grab = subprocess.Popen(['zip', '-sf', file], stdout=subprocess.PIPE)
    grab = grab.communicate()[0]
    grab = grab.decode('-utf8')
    return grab

def get_subfolder(file):
    contents = view_archive(file)
    subfolder = contents.split()[2].strip('\/')
    return subfolder

def unzip_archive(file):
    out = f"unzip {file}"
    print(out)
    subprocess.call(out, shell=True)

def gzip_file(file):
    out = f"gzip {file}"
    print(out)
    subprocess.call(out, shell=True)j

def rm(file):
    out = f"rm {file}"
    print(out)
    subprocess.call(out, shell=True)

if __name__ == '__main__':
    # Kim et al. 10X files
    files = glob.glob('10x*.zip')
    samples = pd.DataFrame({
      'file': files,
      'subfolder': [get_subfolder(f) for f in files]
    })
    samples.to_csv("samples.csv", index=False)
    # unzip each archive
    for f in files:
        unzip_archive(f)
    # gzip each file in a subfolder to save space
    files = glob.glob("*/*")
    for f in files:
        gzip_file(f)
    # SMART-seq data
    unzip_archive("SMART-seq_VMH_cpm.rda.zip")
    files = glob.glob("*.zip")
    for f in files:
        rm(f)

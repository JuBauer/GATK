__author__ = 'jb393'
import os
import sys
import re
import errno
import glob
import time
import pickle
import logging
from joblib import Parallel, delayed
import multiprocessing


def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raisec

def get_filepaths(directory):
    """
    This function will generate the file names in a directory
    tree by walking the tree either top-down or bottom-up. For each
    directory in the tree rooted at directory top (including top itself),
    it yields a 3-tuple (dirpath, dirnames, filenames).
    """
    file_paths = []  # List which will store all of the full filepaths.

    # Walk the tree.
    for root, directories, files in os.walk(directory):
        for filename in files:
            # Join the two strings in order to form the full filepath.
            filepath = os.path.join(root, filename)
            file_paths.append(filepath)  # Add it to the list.

    return file_paths

logFilename = sys.argv[0].split(".")[0]
logging.basicConfig(format='%(asctime)-5s %(message)-10s ',
                    filename=logFilename + ".log",
                    filemode='w',
                    level=logging.DEBUG,
                    datefmt='%m-%d-%Y  %H:%M:%S')

logging.info("*************ChipSeq PIPELINE*************")
logging.info(" ")
logging.info("User command: " + str(sys.argv))


logging.info(" ")
logging.info(" ")
logging.info("#################################")
logging.info("READING PARAMETERS FILE")
params_file = sys.argv[1]
args = {}
with open(params_file, 'r') as f:
    for line in f:
        entry = line.strip().split("=")
        if entry[0]:
            args[entry[0].strip(" ")] = entry[1].strip()

path = args['Working directory'].replace("\\", "")
refGenome = args['Reference Genome']#
nsamples = int(args['Number of samples'])

logging.info("\n -- Input parameters: \n Working Directory = " + path + "\nReference Genome = " + refGenome + "\n Number of samples = " + str(nsamples))
os.chdir(path)
allFiles = get_filepaths(path)
allFiles= [allFiles[i] for i, x in enumerate(allFiles) if re.findall("_L\d{3}.R\d_\d{3}.fastq", x)]

sampleNames = []
logging.info('... Sample names:')
out = open('sample_names.txt', 'w')
for f in allFiles:
    f = f.split('/')[-1]
    f = re.compile("_L\d{3}.R\d_\d{3}.fastq").split(f)[0]
    if f not in sampleNames:
        sampleNames.append(f)
        out.write(f+"\n")
    logging.info(f)
out.close()

# Create rawReads folder and each sample's folder and move reads to corresponding directories
folders = os.listdir('.')
readsFiles = [folders[i] for i, x in enumerate(folders) if re.findall('rawReads',x)]
if not readsFiles:
    logging.info(" ")
    logging.info(" ")
    logging.info("#################################")
    logging.info('... Creating folder structure and organizing reads:')
    make_sure_path_exists('rawReads')
    sampleDir = []
    for sample in sampleNames:
        reads = [allFiles[i] for i,x in enumerate(allFiles) if re.findall(sample,x)]
        if sample not in sampleDir:
            make_sure_path_exists('rawReads/'+sample)
            for r in reads:
               os.system('mv ' + '"' + r + '"' + ' rawReads/' + sample)
               logging.info('Moved reads from sample ' + sample + ' ' + r)
            sampleDir.append(sample)

exit()





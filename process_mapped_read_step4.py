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

def sam2bam(i):
    os.system("samtools view -b alignedReads/" + i + ".sam >alignedReads/"+i+".Aligned.bam")
    logging.info('sam2bam sample'+i+'done')

def indexing(i):
    print("indexing")
    os.system("samtools index alignedReads/dedup/" + i+".Aligned.sortedByCoord.dedup.out.bam" )
    logging.info('Sample ' + i + ' done')

def bamstats(i):
    os.system("~/bin/bamtools stats -in alignedReads/" + i + ".Aligned.bam >bamstats/"+i+".bam_stats.txt")
    logging.info('Sample ' + i + ' bamtools stats, done')

def sortByCoord(i):
    os.system('samtools sort alignedReads/' + i + '.Aligned.bam alignedReads/' + i + '.Aligned.sortedByCoord.out')
    logging.info("Sample " + i + ' done')

def removeDupli(i):
    os.system('java -jar ~/bin/picard.jar MarkDuplicates INPUT=alignedReads/' + i + '.Aligned.sortedByCoord.out.bam OUTPUT=alignedReads/dedup/' + i + '.Aligned.sortedByCoord.dedup.out.bam REMOVE_DUPLICATES=F METRICS_FILE=alignedReads/dedup/'+i+'_dedup_metrics.txt')
    logging.info("Sample " + i + ' done')

logFilename = './' + sys.argv[0].split(".")[0].split('/')[-1]
logging.basicConfig(format='%(asctime)-5s %(message)-10s ',
                    filename=logFilename + ".log",
                    #filemode='w',
                    level=logging.DEBUG,
                    datefmt='%m-%d-%Y  %H:%M:%S')

logging.info(" ")
logging.info(" ")
logging.info("***************************************")
logging.info("*************Processing Mapped Reads*************")
logging.info(" ")
logging.info("User command: " + str(sys.argv))


logging.info(" ")
logging.info(" ")
logging.info("#################################")
logging.info("... READING PARAMETERS FILE")
params_file = sys.argv[1]
args = {}
with open(params_file, 'r') as f:
    for line in f:
        entry = line.strip().split("=")
        if entry[0]:
            args[entry[0].strip(" ")] = entry[1].strip()
logging.info('-- Input parameters:')
path = args['Working directory'].replace("\\", "")
logging.info('Working Directory = ' + path)

nsamples = int(args['Number of samples'])
os.chdir(path)

logging.info(" ")
logging.info(" ")
logging.info("#################################")
# Read in sampleNames
logging.info('... Sample names:')
sampleNames = []
sample_names_file = open('sample_names.txt','r')
for line in sample_names_file:
    sampleNames.append(line.strip())
    logging.info(line.strip())


# Index Bam files
logging.info(" ")
logging.info(" ")
logging.info("#################################")
logging.info("... SAM 2 BAM files")
alignmentOut = os.listdir(path+"/alignedReads/") # Check if aligned reads files exist
bamFiles = get_filepaths(path+"/alignedReads/")
bamFiles = [bamFiles[i] for i, x in enumerate(alignmentOut) if re.findall("\.bam", x)] # Check if bam files exist
sampleToBam=list()
for sample in sampleNames:
    found=False
    for files in bamFiles:
        file=files.split(".")[0]
        file2=file.split("/")[-1]
      #  print(file2+"<->"+sample)
        if sample == file2:
            found=True
    if not found:
        logging.info(sample)
        sampleToBam.append(sample)
        print(sample)

#sys.exit("Done")
if len(sampleToBam)>0:
    Parallel(n_jobs=8)(delayed(sam2bam)(i) for i in sampleToBam)
logging.info("... Finished SAM 2 BAM")


# Mapping QC
logging.info(" ")
logging.info(" ")
logging.info("#################################")
logging.info("... Mapping QC")
qcOut = os.listdir(path+"/alignedReads/")
make_sure_path_exists(path+"/bamstats/")
bamstatsOut=os.listdir(path+"/bamstats/")
indicesQCFiles=get_filepaths(path+"/bamstats")
indicesQCFiles = [indicesQCFiles[i] for i, x in enumerate(bamstatsOut) if re.findall(".bam_stats.txt", x)] # Check if aligned reads files exist
sampleToBamStat=list()
for sample in sampleNames:
    found=False
    for files in indicesQCFiles:
        #print(files)
        file=files.split(".")[0]
        file2=file.split("/")[-1]
        if sample == file2:
            found=True
    if not found:
        logging.info(sample)
        sampleToBamStat.append(sample)

if len(sampleToBamStat)>0:
    Parallel(n_jobs=8)(delayed(bamstats)(i) for i in sampleToBamStat)
logging.info("... Finished generating mapping QC")
logging.info('Mapping summary done')
logging.info(" ")

#sort by coordinate
logging.info(" ")
logging.info(" ")
logging.info("#################################")
logging.info("... Sorting BAM files by coord ")
# Sort bam files by name
alignmentOut = os.listdir(path+"/alignedReads/") # Check if aligned reads files exist
sortedFile=get_filepaths(path+"/alignedReads")
sortedFile = [sortedFile[i] for i, x in enumerate(alignmentOut) if re.findall(".sortedByCoord.out.bam", x)]
bamToSort=list()
for sample in sampleNames:
    found=False
    for files in sortedFile:
        #print(files)
        file=files.split(".")[0]
        file2=file.split("/")[-1]
        if sample == file2:
            found=True
    if not found:
        logging.info(sample)
        bamToSort.append(sample)

if len(bamToSort)>0:
    Parallel(n_jobs=8)(delayed(sortByCoord)(i) for i in bamToSort)
logging.info("... Finished sorting BAM files by name")


# Removing duplicate
logging.info(" ")
logging.info(" ")
logging.info("#################################")
logging.info("... remove duplicate")
make_sure_path_exists(path+"/alignedReads/dedup/")
dedupOut = os.listdir(path+"/alignedReads/")
indicesdedupFiles = [i for i, x in enumerate(dedupOut) if re.findall(".dedup*", x)] # Check if aligned reads files exist
if not indicesdedupFiles:
    Parallel(n_jobs=8)(delayed(removeDupli)(i) for i in sampleNames)
logging.info("... Finished removing duplicates")
logging.info('Mapping summary done')
logging.info(" ")

# Index Bam files
logging.info(" ")
logging.info(" ")
logging.info("#################################")
logging.info("... Indexing BAM files")
indexedOut = os.listdir(path+"/alignedReads/")
indexedFile=get_filepaths(path+"/alignedReads")
indexedFile = [indexedFile[i] for i, x in enumerate(indexedOut) if re.findall(".bai", x)]
bamToIndex=list()
for sample in sampleNames:
    found=False
    for files in indexedFile:
        #print(files)
        file=files.split(".")[0]
        file2=file.split("/")[-1]
        if sample == file2:
            found=True
    if not found:
        logging.info(sample)
        bamToIndex.append(sample)

if len(bamToIndex)>0:
    Parallel(n_jobs=8)(delayed(indexing)(i) for i in bamToIndex)
logging.info("... Finished Indexing")

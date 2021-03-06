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

#################### ~~~ Realign function ~~~~ #########
def realign(i):
    os.system("java -Xmx1G -jar ~/GenomeAnalysisTK.jar --consensusDeterminationModel KNOWNS_ONLY -LOD 0.4 -T IndelRealigner -R "+refGenome+" -targetIntervals "+path+"/realigned/"+refGenome+".intervals -i "+path+"alignedReads/dedup/"+i+".Aligned.sortedByCoord.dedup.out.bam -o "+path+"/realigned/"+i+".realignedBam.bam  -known "+vcfIndel+" >realigned/"+i+".log")
    logging.info('target realignement on sample '+i+' done')


#################### ~~~ Base recalibration function ~~~~ #########
def recal(i):
    os.system("java -Xmx1G -jar ~/GenomeAnalysisTK.jar -T PrintReads -R "+refGenome+" -i "+path+"/realigned/"+i+".realignedBam.bam -o "+path+"/recalibration/"+i+".recal.bam  -BQSR "+path+"recalibration/"+i+"_report.grp")
    logging.info('target realignement on sample '+i+' done')

###################### ~~~~ INIT log file ~~~~~ ####################
logFilename = './' + sys.argv[0].split(".")[0].split('/')[-1]
logging.basicConfig(format='%(asctime)-5s %(message)-10s ',
                    filename=logFilename + ".log",
                    #filemode='w',
                    level=logging.DEBUG,
                    datefmt='%m-%d-%Y  %H:%M:%S')

####################################################################

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


################ ~~~~ READING PARAMETERS ~~~~ ##################
##### NEEDED : Working Directory, Reference Genome, VCF indels, VCF snps

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
refGenome = args['Reference Genome']
logging.info('Reference Genome = '+refGenome)
vcfIndel = args['VCF indels']
vcfSNP = args['VCF snps']
logging.info('VCF indels ref = '+vcfIndel)
logging.info('VCF snps ref = '+vcfSNP)

nsamples = int(args['Number of samples'])
os.chdir(path)

logging.info(" ")
logging.info(" ")
logging.info("#################################")

########## ~~~~~ Read in Samples Names ~~~~~ ###############
#### NEEDED sample_names.txt created in step 2

logging.info('... Sample names:')
sampleNames = []
sample_names_file = open('sample_names.txt','r')
for line in sample_names_file:
    sampleNames.append(line.strip())
    logging.info(line.strip())


############# ~~~~~ START Realigner ~~~~~ ################
logging.info(" ")
logging.info(" ")
logging.info("#################################")
logging.info("... Indel realigner, create target file")

############ ~~~~ create target file ~~~~~ ###############
##create output folder (check if it exists)
make_sure_path_exists(path+"/realigned/")
##check if target file already exists
realignmentOut = os.listdir(path+"/realigned/")
intervalFiles = get_filepaths(path+"/realigned/")
intervalFiles = [intervalFiles[i] for i, x in enumerate(realignmentOut) if re.findall("\.intervals", x)]
#sys.exit("Done")
if len(intervalFiles)>0:
    os.system("java -Xmx1G -jar ~/GenomeAnalysisTK.jar -T RealignerTargetCreator -R "+refGenome+" -o "+path+"/realigned/"+refGenome+".intervals  -known "+vcfIndel+" >realigned/target_creation.log")
logging.info("... Finished interval target creation")




####################~~~~  Realign individual sample~~~~~ ########
logging.info(" ")
logging.info(" ")
logging.info("#################################")
logging.info("... Realignement for each sample")
make_sure_path_exists(path+"/recalibration/")
realignOut = os.listdir(path+"/recalibration/")
realignFiles=get_filepaths(path+"/recalibration")
realignFiles = [realignFiles[i] for i, x in enumerate(realignOut) if re.findall("realignedBam.bam", x)] # Check if realigned reads files exist
sampleToRealign=list()
for sample in sampleNames:
    found=False
    for files in realignFiles:
        file=files.split(".")[0]
        file2=file.split("/")[-1]
        if sample == file2:
            found=True
    if not found:
        logging.info(sample)
        sampleToRealign.append(sample)

if len(sampleToRealign)>0:
    Parallel(n_jobs=8)(delayed(recal)(i) for i in sampleToRealign)
logging.info("... Finished sample realignment")
logging.info(" ")

################ ~~~~ Base Recalibration ~~~~ #########################
logging.info(" ")
logging.info(" ")
logging.info("#################################")
logging.info("... Base recalibration")
# Sort bam files by name
recalOut = os.listdir(path+"/recalibration/")
recalFiles=get_filepaths(path+"/recalibration")
recalFiles = [recalFiles[i] for i, x in enumerate(recalOut) if re.findall("recalBam.bam", x)] # Check if realigned reads files exist
sampleToRecal=list()
for sample in sampleNames:
    found=False
    for files in recalFiles:
        file=files.split(".")[0]
        file2=file.split("/")[-1]
        if sample == file2:
            found=True
    if not found:
        logging.info(sample)
        sampleToRecal.append(sample)

if len(sampleToRecal)>0:
    Parallel(n_jobs=8)(delayed(recal)(i) for i in sampleToRecal)
logging.info("... Finished sample base recalibration")
logging.info(" ")



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

#################### ~~~ Call variant function ~~~~ #########
def callVariant(i):
    if len(ploidyDict)>0:
        os.system("java -Xmx1G -jar ~/GenomeAnalysisTK.jar -T UnifiedGenotyper -R "+refGenome+" -i "+path+"recalibration/"+i+".recal.bam -o "+path+"/variant/"+i+".raw.vcf  -dbsnp "+vcfSNP+" -stand_call_conf [50.0] -stand_emit_conf 10.0 -L "+bedfile+" -ploidy "+ploidyDict[i]+" >variant/"+i+"raw.log")
    else:
        os.system("java -Xmx1G -jar ~/GenomeAnalysisTK.jar -T UnifiedGenotyper -R "+refGenome+" -i "+path+"recalibration/"+i+".recal.bam -o "+path+"/variant/"+i+".raw.vcf  -dbsnp "+vcfSNP+" -stand_call_conf [50.0] -stand_emit_conf 10.0 -L "+bedfile+" >variant/"+i+"raw.log")
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
##### NEEDED : Working Directory, Reference Genome, vcfSNP, BED File, Variable ploidy, Ploidy Table (if variable ploidy TRUE)

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
vcfSNP = args['VCF snps']
bedFile=args['BED File']
ploidy= args['Variable ploidy']
ploidyDict={}
if not ploidy=="FALSE":
    ploidyFile= args['Ploidy Table']
    for ploidy in ploidyFile:
        ploidy = ploidy.strip()
        ploidyTab = ploidy.split("\t")
        ploidyDict[ploidyTab[0]]=ploidyTab[1]

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

############# ~~~~~ START variant calling ~~~~~ ################
logging.info(" ")
logging.info(" ")
logging.info("#################################")
logging.info("... Variant calling")

############ ~~~~ create target file ~~~~~ ###############
##create output folder (check if it exists)
make_sure_path_exists(path+"/variant/")
##check if target file already exists
variantOut = os.listdir(path+"/variant/")
variantFiles = get_filepaths(path+"/variant/")
variantFiles = [variantFiles[i] for i, x in enumerate(variantOut) if re.findall("\.vcf", x)]
sampleToCall=list()
for sample in sampleNames:
    found=False
    for files in variantFiles:
        file=files.split(".")[0]
        file2=file.split("/")[-1]
        if sample == file2:
            found=True
    if not found:
        logging.info(sample)
        sampleToCall.append(sample)

if len(sampleToCall)>0:
    Parallel(n_jobs=8)(delayed(callVariant)(i) for i in sampleToCall)
logging.info("... Finished variant calling ")
logging.info(" ")
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

def alignment(i):
    print("aligning2?")
    trimmedReads = os.listdir("./readyToMap")
    trimmedReads = [trimmedReads[y] for y, x in enumerate(trimmedReads) if re.findall(i, x)]
    r2 = [trimmedReads[y] for y, x in enumerate(trimmedReads) if re.findall("_all_lane_R2.fastq", x)]
    if r2:
        r1 = [trimmedReads[y] for y, x in enumerate(trimmedReads) if re.findall( "_all_lane_R1.fastq", x)]
        r2 = [trimmedReads[y] for y, x in enumerate(trimmedReads) if re.findall("_all_lane_R2.fastq", x)]
        print(r1[0]+" "+r2[0])
        os.system("bwa mem  -t 8 "+ refGenome +" "+path+"/readyToMap/" + r1[0] + " "+path+"/readyToMap/" + r2[0] + " >"+path+"/alignedReads/" + i+".sam")
        logging.info('Sample ' + i + ' done, pairedEnd mode')
    else:
        r1 = [trimmedReads[y] for y, x in enumerate(trimmedReads) if re.findall( "_all_lane_R1.fastq", x)]
        os.system("bwa mem  -t 8 "+ refGenome +" "+path+"/readyToMap/" + r1[0] + ' >'+path+'/alignedReads/' + i+".sam")
        logging.info('Sample ' + i + ' done, singleEnd mode')

def sam2bam(i):
    os.system("samtools view -b alignedReads/" + i + ".sam >alignedReads/"+i+".Aligned.bam")
    logging.info('sam2bam sample'+i+'done')

def indexing(i):
    print("indexing")
    os.system("samtools index alignedReads/1pass_ann/" + i+"Aligned.sortedByCoord.dedup.out.bam" )
    logging.info('Sample ' + i + ' done')

def bamstats(i):
    os.system("~/bin/bamtools stats -in alignedReads/" + i + ".Aligned.bam >bamstats/"+i+"_bam_stats.txt")
    logging.info('Sample ' + i + ' bamtools stats, done')

def bam2bed(i):
    os.system("bedtools bamtobed -i alignedReads/1pass_ann/" + i + "Aligned.sortedByCoord.out.bam >alignedReads/1pass_ann/"+i+".bed")
    logging.info('Sample ' + i + ' bedtools bamToBed, done')

def junctions(i):
    os.system("python ~/bin/junction_annotation.py -i alignedReads/1pass_ann/" + i +"Aligned.sortedByCoord.out.bam  -o alignedReads/1pass_ann/QC/" + i + " -r " + bedFile + " & ")
    os.system("python ~/bin/junction_saturation.py -i alignedReads/1pass_ann/" + i +"Aligned.sortedByCoord.out.bam  -o alignedReads/1pass_ann/QC/" + i + " -r " + bedFile)
    logging.info("Calculated junction annotation and junction saturation sample " + i)

def sortByName(i):
    os.system('samtools sort -n alignedReads/1pass_ann/' + i + 'Aligned.sortedByCoord.out.bam alignedReads/1pass_ann/' + i + 'Aligned.sortedByCoord.sortedByName.out')
    logging.info("Sample " + i + ' done')

def sortByCoord(i):
    os.system('samtools sort alignedReads/1pass_ann/' + i + 'Aligned.bam alignedReads/1pass_ann/' + i + 'Aligned.sortedByCoord.out')
    logging.info("Sample " + i + ' done')

def removeDupli(i):
    os.system('java -jar ~/bin/picard.jar MarkDuplicates INPUT=alignedReads/1pass_ann/' + i + 'Aligned.sortedByCoord.out.bam OUTPUT=alignedReads/1pass_ann/' + i + 'Aligned.sortedByCoord.dedup.out.bam REMOVE_DUPLICATES=T METRICS_FILE=alignedReads/1pass_ann/'+i+'_dedup_metrics.txt')
    logging.info("Sample " + i + ' done')


###################
#  * Rscripts used:
# mapping_summary.R
# mapping_distribution.R
###################

logFilename = './' + sys.argv[0].split(".")[0].split('/')[-1]
logging.basicConfig(format='%(asctime)-5s %(message)-10s ',
                    filename=logFilename + ".log",
                    #filemode='w',
                    level=logging.DEBUG,
                    datefmt='%m-%d-%Y  %H:%M:%S')

logging.info(" ")
logging.info(" ")
logging.info("***************************************")
logging.info("*************MAPPING READS*************")
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
refGenome = args['Reference Genome']
logging.info('Reference Genome = ' + refGenome)
nsamples = int(args['Number of samples'])
#logging.info("\n -- Input parameters: \n Working Directory = " + path + "\n GTF file = " + gtfFile + "\n Reference Genome = " + refGenome + "\n Number of samples = " + str(nsamples) + "\n BedFile = " + bedFile + "\n BedFile_10k = " + bedFile_10k + "\n refFlat = " + refFlat + "\n rRNA_interval_list = " + rRNA_interval_list + "\n strand = " + strand )
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


logging.info(" ")
logging.info(" ")
logging.info("#################################")
logging.info("... Mapping of reads")
make_sure_path_exists("./alignedReads/")
alignmentOut = os.listdir(path+"/alignedReads/") # Check if aligned reads files exist
indicesAlignmentFiles = [i for i, x in enumerate(alignmentOut) if re.findall(".sam", x)] # Check if aligned reads files exist

if not indicesAlignmentFiles:
    print("aligning?")
    Parallel(n_jobs=2)(delayed(alignment)(i) for i in sampleNames)
logging.info("... Finished Mapping")




# logging.info(" ")
# logging.info(" ")
# logging.info("#################################")
# logging.info("... Convert bam to bed")
# bedOut = os.listdir(path+"/alignedReads/1pass_ann")
# indicesbedFiles = [i for i, x in enumerate(bedOut) if re.findall(".bed", x)] # Check if aligned reads files exist
# if not indicesbedFiles:
#     Parallel(n_jobs=8)(delayed(bam2bed)(i) for i in sampleNames)
# logging.info("... Finished converting bam to bed")
# logging.info('Mapping summary done')
# logging.info(" ")



logging.info(" ")
logging.info(" ")
logging.info(" ")
logging.info("##################################################################")
logging.info("------------------------------------------------------------------")
logging.info("##################################################################")
logging.info(" ")
logging.info(" ")
logging.info(" ")
# GeneBody Coverage
# make_sure_path_exists("./alignedReads/1pass_ann/QC")
# make_sure_path_exists("./Report/figure/")
# alignmentQCOut = os.listdir(path+"/alignedReads/1pass_ann/QC") # Check if aligned reads files exist
# indicesAlignmentQCFiles = [i for i, x in enumerate(alignmentQCOut) if re.findall(".geneBodyCoverage.", x)] # Check if aligned reads files exist
# if not indicesAlignmentQCFiles:
#     os.system("ls ./alignedReads/1pass_ann/*Aligned.sortedByCoord.out.bam > tempbamfiles.txt")
#     os.system("python ~/bin/geneBody_coverage.py -r " + bedFile_10k + " -i tempbamfiles.txt -o alignedReads/1pass_ann/QC/10KGenes")
#     os.system("rm tempbamfiles.txt")
# os.system('cp alignedReads/1pass_ann/QC/10KGenes.geneBodyCoverage.curves.pdf Report/figure/10KGenes_geneBodyCoverage_curves.pdf')
# logging.info('GeneBody Coverage done')
# logging.info(" ")
#
# # Junctions and junction saturation
# qcOut = os.listdir(path+"/alignedReads/1pass_ann/QC")
# qcFiles = [i for i, x in enumerate(qcOut) if re.findall(".junction.", x)] # Check if aligned reads files exist
# if not qcFiles:
#      Parallel(n_jobs=8)(delayed(junctions)(i) for i in sampleNames)
# os.system('grep "y=c(" alignedReads/1pass_ann/QC/*junctionSaturation*  | sed \'s/:y=c(/,/g\' | sed \'s/.junctionSaturation_plot.r//g\' | sed \'s/)//g\' | sed \"s/.*\///g\"  > alignedReads/1pass_ann/QC/junctionSat_all.csv')
# os.system('~/bin/junctionPlotAll.R .')
# logging.info(" ")
#
# # Collect Metrics
# qcOut = os.listdir(path+"/alignedReads/1pass_ann/QC") # Check if aligned reads files exist
# qcFiles = [i for i, x in enumerate(qcOut) if re.findall(".metrics.txt", x)] # Check if aligned reads files exist
# logging.info("Calculating RNASeq Metrics")
# if not qcFiles:
#     os.system("find alignedReads/1pass_ann/*.sortedByCoord.out.bam | sed \"s/.*\///g\" | sed \"s/.sortedByCoord.out.bam//g\" | " # get file names and format them
#           "parallel -j 3 --no-notice "
#           "\"java -jar ~/tools/picard-tools-1.127/picard.jar CollectRnaSeqMetrics "
#           "REF_FLAT=" + refFlat + " "
#           "RIBOSOMAL_INTERVALS=" + rRNA_interval_list + " "
#           "STRAND_SPECIFICITY=" + strand + " "
#           "INPUT=alignedReads/1pass_ann/{}.sortedByCoord.out.bam "
#           "OUTPUT=alignedReads/1pass_ann/QC/{}_metrics.txt \"")
#     logging.info("Calculated RNASeq Metrics")
# def pct(i):
#     os.system('mapping_distribution.R ./alignedReads/1pass_ann/QC/ ' + i)
# Parallel(n_jobs=8)(delayed(pct)(i) for i in sampleNames)
# logging.info("... Finished mapping QC")
#




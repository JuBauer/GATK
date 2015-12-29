#!usr/bin/python
#!usr/local/Rscript
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

#if need to rerun preprocessing, remove created file in the rawReads folder
def preprocessing(i):
	allFiles = os.listdir("./rawReads/" + i )
	pairedReads_temp = [allFiles[y] for y, x in enumerate(allFiles) if re.findall("_L\d{3}.R2", x)]
	os.system("mkdir rawReads/" + i + "/tmp")
	os.system("cat ./rawReads/" + i + "/*R1*" + " > ./rawReads/" + i + "/tmp/" + i + "_all_lane_R1.fastq" + gz)
	os.system("mv ./rawReads/" + i + "/tmp/*_all_lane_R1.fastq"+ gz + " ./rawReads/" + i + " && rm -r ./rawReads/" + i + "/tmp")
	os.system("fastqc ./rawReads/" + i + "/" + i + "_all_lane_R1.fastq" + gz + " --outdir=rawReads/" + i+' &>log_fastqc_'+i+'.txt')
	os.system("unzip ./rawReads/" + i + "/" + i + "_all_lane_R1_fastqc.zip -d rawReads/" + i + "/")
	if pairedReads_temp:
		os.system("mkdir rawReads/" + i + "/tmp")
		os.system("cat ./rawReads/" + i + "/*R2*" + " > ./rawReads/" + i + "/tmp/" + i + "_all_lane_R2.fastq" + gz)
		os.system("mv ./rawReads/" + i + "/tmp/*_all_lane_R2.fastq"+ gz + " ./rawReads/" + i + " && rm -r rawReads/" + i + "/tmp")
		os.system("fastqc ./rawReads/" + i + "/" + i + "_all_lane_R2.fastq" + gz + " --outdir=rawReads/" + i+' &>log_fastqc_'+i+'.txt')
		os.system("unzip ./rawReads/" + i + "/" + i + "_all_lane_R2_fastqc.zip -d rawReads/" + i + "/")
		os.chdir(path)

	#logging.info("Preprocessing done for sample " + i)

def tables(i):
	os.system('python ~/bin/fastqc_plot_data.py '+path+'/rawReads/' + i + '/' + i + '_all_lane_R1_fastqc/fastqc_data.txt all ./rawReads/' + i + '/' + i + '_all_lane_R1_fastqc/')
	allFiles = os.listdir("./rawReads/" + i )
	pairedReads_temp = [allFiles[y] for y, x in enumerate(allFiles) if re.findall("_L\d{3}.R2", x)]
	if pairedReads_temp:
		os.system('python ~/bin/fastqc_plot_data.py '+path+'/rawReads/' + i + '/' + i + '_all_lane_R2_fastqc/fastqc_data.txt all '+path+'/rawReads/' + i + '/' + i + '_all_lane_R2_fastqc/')
		os.system('python ~/bin/fastqc_plot_data.py '+path+'/trimmedReads/' + i + '_all_lane_R1_val_1_fastqc/fastqc_data.txt all '+path+'/trimmedReads/' + i + '_all_lane_R1_val_1_fastqc/')
		os.system('python ~/bin/fastqc_plot_data.py '+path+'/trimmedReads/' + i + '_all_lane_R2_val_2_fastqc/fastqc_data.txt all '+path+'/trimmedReads/' + i + '_all_lane_R2_val_2_fastqc/')
	else:
		os.system('python ~/bin/fastqc_plot_data.py '+path+'/trimmedReads/' + i + '_all_lane_R1_trimmed_fastqc/fastqc_data.txt all '+path+'/trimmedReads/' + i + '_all_lane_R1_trimmed_fastqc/')

def moveFiles(i):
	allFiles=os.listdir("./rawReads/"+i)
	mergedReads=[allFiles[y] for y,x in enumerate(allFiles) if re.findall("all_lane.+.gz", x)]
	if mergedReads:
		os.system("mv "+path+"/rawReads/"+i+"/"+i+"_all_lane_* "+path+"/readyToMap/")
		os.system("gunzip "+path+"/readyToMap/*.gz")




###################
#  * Rscripts used:
# preprocessing_numbers.R
# trimgalore_summary.R
# mapping_summary.R
###################

logFilename = './' + sys.argv[0].split(".")[0].split('/')[-1]
logging.basicConfig(format='%(asctime)-5s %(message)-10s ',
					filename=logFilename + ".log",
					#filemode='w',
					level=logging.DEBUG,
					datefmt='%m-%d-%Y  %H:%M:%S')


logging.info("*************PREPROCESSING-READS*************")
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
logging.info('Number of samples = ' + str(nsamples))
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

# Check if cat files exist, might be removed
catFiles = []
for root, dir, files in os.walk(path):
	catFiles.extend(files)
indicesCatFiles = [catFiles[i] for i, x in enumerate(catFiles) if re.findall("all_lane", x)]

# Check if input files are uncompressed or .gz
readFiles = []
for root, dir, files in os.walk(path):
	readFiles.extend(files)
indicesgzFiles = [i for i, x in enumerate(readFiles) if re.findall(".fastq.gz", x)]
logging.info(' ')
if indicesgzFiles:
   gz =".gz"
   logging.info("Input files are .gz ")
else:
   gz =""
   logging.info("Input files are uncompressed ")

logging.info(" ")
logging.info(" ")
logging.info("#################################")
logging.info("... PRE-PROCESSING OF READS")
make_sure_path_exists("./trimmedReads")

# Running jobs
if not indicesCatFiles:
	Parallel(n_jobs=8)(delayed(preprocessing)(i) for i in sampleNames)


logging.info(" ")
logging.info(" ")
logging.info("#################################")
logging.info("... MOVING READS")
make_sure_path_exists(path+"/readyToMap")
readyReadFiles=[]
for root, dir, files in os.walk(path+"/readyToMap"):
	readyReadFiles.extend(files)
readyReadFiles = [readyReadFiles[i] for i, x in enumerate(readyReadFiles) if re.findall("all_lane", x)]

if not readyReadFiles:
	Parallel(n_jobs=8)(delayed(moveFiles)(i) for i in sampleNames)


logging.info(" ")
logging.info(" ")
logging.info("#################################")
logging.info("... Creating files with output information and summary")
# Create files with output information (used in the Report)
os.chdir(path)
#os.system("~/bin/preprocessing_numbers.R .")
#os.system("~/bin/trimgalore_summary.R .")

# Create files from fastqc_data.txt for creating plots (used in the Report)

Parallel(n_jobs=8)(delayed(tables)(i) for i in sampleNames)

os.system('ls rawReads/*/*fastqc  |  grep -v trimmed  | grep ":"  | sed \'s/://g\' > sample_names2.txt')
#os.system('python ~/bin/fastqc_summary.py ./sample_names2.txt ./summary_fastqc.txt')

logging.info(" ")
logging.info(" ")
logging.info(" ")
logging.info("##################################################################")
logging.info("------------------------------------------------------------------")
logging.info("##################################################################")
logging.info(" ")
logging.info(" ")
logging.info(" ")
#logging.info("DONE PRE-PROCESSING OF READS")


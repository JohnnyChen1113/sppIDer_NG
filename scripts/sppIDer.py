__author__ = 'Quinn'
__modified_by__ = 'Junhao Chen'

import argparse, multiprocessing, sys, re, subprocess, time
import os

################################################################
# This script runs the full sppIDer pipeline.
# Before running this you must make your combination reference genome with the script combineRefGenomes.py
# This script will map short-read data to a combination reference genome and parse the outputs to create a summary of
# where and how well the reads map to each species in the combination reference genome.
#
# Program Requirements: bwa, samtools, bedtools, R, Rpackage ggplot2, Rpackage dplyr
# Input: Output name, Combination reference genome, fastq short-read files
#
################################################################

parser = argparse.ArgumentParser(description="Run full sppIDer")
parser.add_argument('--out', help="Output prefix, required", required=True)
parser.add_argument('--ref', help="Reference Genome, required", required=True)
parser.add_argument('--r1', help="Read1, required", required=True)
parser.add_argument('--r2', help="Read2, optional")
parser.add_argument('--byBP', help="Calculate coverage by basepair, optional, DEFAULT, can't be used with --byGroup", dest='bed', action='store_true')
parser.add_argument('--byGroup', help="Calculate coverage by chunks of same coverage, optional, can't be used with --byBP", dest='bed', action='store_false')
parser.add_argument('--seq-name', help="Output sequence name, optional")
parser.add_argument('--mq-threshold', type=int, help="Set MQ threshold, optional", default=3)
parser.add_argument('--cores', type=int, help="Set number of cores used in analysis, optional", default=(multiprocessing.cpu_count()/2))
parser.add_argument('--seq-type', choices=['PacBio', 'ONT'], help="Set sequence type (PacBio, ONT), optional")
parser.add_argument('--mapping-tool', choices=['bwa', 'minimap2'], help="Set mapping tool (default is bwa, optional is minimap2)", default='bwa')
parser.set_defaults(bed=True)
args = parser.parse_args()

# docker vars
scriptDir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '')
workingDir = os.path.join(os.getcwd(), '')

numCores = str(args.cores)
outputPrefix = args.out
refGen = args.ref
read1Name = args.r1
read2Name = args.r2 if args.r2 else None
start = time.time()

def calcElapsedTime(endTime):
    trackedTime = str()
    if 60 < endTime < 3600:
        min = int(endTime) // 60
        sec = int(endTime % 60)
        trackedTime = f"{min} mins {sec} secs"
    elif 3600 < endTime < 86400:
        hr = int(endTime) // 3600
        min = int((endTime % 3600) // 60)
        sec = int(endTime % 60)
        trackedTime = f"{hr} hrs {min} mins {sec} secs"
    elif 86400 < endTime < 604800:
        day = int(endTime) // 86400
        hr = int((endTime % 86400) // 3600)
        min = int((endTime % 3600) // 60)
        sec = int(endTime % 60)
        trackedTime = f"{day} days {hr} hrs {min} mins {sec} secs"
    elif 604800 < endTime:
        week = int(endTime) // 604800
        day = int((endTime % 604800) // 86400)
        hr = int((endTime % 86400) // 3600)
        min = int((endTime % 3600) // 60)
        sec = int(endTime % 60)
        trackedTime = f"{week} weeks {day} days {hr} hrs {min} mins {sec} secs"
    else:
        trackedTime = f"{int(endTime)} secs"
    return trackedTime

trackerOut = open(os.path.join(workingDir, outputPrefix + "_sppIDerRun.info"), 'w')
trackerOut.write(f"outputPrefix={args.out}\n")
trackerOut.write(f"ref={refGen}\n")
trackerOut.write(f"read1={read1Name}\n")
if read2Name: trackerOut.write(f"read2={read2Name}\n")
if not args.bed:
    trackerOut.write("coverage analysis option = by coverage groups, bedgraph format -bga\n")
else: trackerOut.write("coverage analysis option = by each base pair -d\n")
trackerOut.close()

########################## BWA ###########################
#bwaOutName = outputPrefix + ".sam"
#bwaOutFile = open(os.path.join(workingDir, bwaOutName), 'w')
#if read2Name:
#    print(f"Read1={read1Name}\nRead2={read2Name}")
#    subprocess.call([args.mapping_tool, "mem", "-t", numCores, refGen, read1Name, read2Name], stdout=bwaOutFile, cwd=workingDir)
#else:
#    print(f"Read1={read1Name}")
#    subprocess.call([args.mapping_tool, "mem", "-t", numCores, refGen, read1Name], stdout=bwaOutFile, cwd=workingDir)
#bwaOutFile.close()
#print("BWA complete")
#currentTime = time.time() - start
#elapsedTime = calcElapsedTime(currentTime)
#print(f"Elapsed time: {elapsedTime}")
#trackerOut = open(os.path.join(workingDir, outputPrefix + "_sppIDerRun.info"), 'a')
#trackerOut.write(f"BWA complete\nElapsed time: {elapsedTime}")
#trackerOut.close()

bwaOutName = outputPrefix + ".sam"
bwaOutFile = open(os.path.join(workingDir, bwaOutName), 'w')

match args.mapping_tool:
    case 'bwa':
        if args.seq_type and args.seq_type.lower() == 'pacbio':
            subprocess.call([args.mapping_tool, "mem", "-t", numCores, "-x", "pacbio", refGen, read1Name], stdout=bwaOutFile, cwd=workingDir)
        elif args.seq_type and args.seq_type.lower() == 'ont':
            subprocess.call([args.mapping_tool, "mem", "-t", numCores, "-x", "ont2d", refGen, read1Name], stdout=bwaOutFile, cwd=workingDir)
        else:
            subprocess.call([args.mapping_tool, "mem", "-t", numCores, refGen, read1Name], stdout=bwaOutFile, cwd=workingDir)
    case 'minimap2':
        if args.seq_type and args.seq_type.lower() == 'pacbio':
            subprocess.call(["minimap2", "-x", "map-pb", "-a", "-t", numCores, refGen, read1Name], stdout=bwaOutFile, cwd=workingDir)
        elif args.seq_type and args.seq_type.lower() == 'ont':
            subprocess.call(["minimap2", "-x", "map-ont", "-a", "-t", numCores, refGen, read1Name], stdout=bwaOutFile, cwd=workingDir)
        else:
            subprocess.call(["minimap2", "-a", "-t", numCores, refGen, read1Name], stdout=bwaOutFile, cwd=workingDir)
    case _:
        raise ValueError(f"Unsupported mapping tool: {args.mapping_tool}")

bwaOutFile.close()
print("BWA complete")
currentTime = time.time() - start
elapsedTime = calcElapsedTime(currentTime)
print(f"Elapsed time: {elapsedTime}")
trackerOut = open(os.path.join(workingDir, outputPrefix + "_sppIDerRun.info"), 'a')
trackerOut.write(f"BWA complete\nElapsed time: {elapsedTime}")
trackerOut.close()

########################## samtools ###########################
samViewOutQual = outputPrefix + ".view.bam"
bamSortOut = outputPrefix + ".sort.bam"
samViewQualFile = open(os.path.join(workingDir, samViewOutQual), 'w')
subprocess.call(["samtools", "view", "-@", numCores, "-q", str(args.mq_threshold), "-bhSu", bwaOutName], stdout=samViewQualFile, cwd=workingDir)
samViewQualFile.close()
subprocess.call(["samtools", "sort", "-@", numCores, samViewOutQual, "-o", bamSortOut], cwd=workingDir)
print("SAMTOOLS complete")
currentTime = time.time() - start
elapsedTime = calcElapsedTime(currentTime)
print(f"Elapsed time: {elapsedTime}")
trackerOut = open(os.path.join(workingDir, outputPrefix + "_sppIDerRun.info"), 'a')
trackerOut.write(f"\nSAMTOOLS complete\nElapsed time: {elapsedTime}")
trackerOut.close()

########################## parse SAM file ###########################
subprocess.call(["python3", os.path.join(scriptDir, "parseSamFile.py"), outputPrefix], cwd=workingDir)
print("Parsed SAM file")
currentTime = time.time() - start
elapsedTime = calcElapsedTime(currentTime)
print(f"Elapsed time: {elapsedTime}")
trackerOut = open(os.path.join(workingDir, outputPrefix + "_sppIDerRun.info"), 'a')
trackerOut.write(f"\nParsed SAM\nElapsed time: {elapsedTime}")
trackerOut.close()

########################## plot MQ scores ###########################
subprocess.call(["Rscript", os.path.join(scriptDir, "MQscores_sumPlot.R"), outputPrefix], cwd=workingDir)
print("Plotted MQ scores")
currentTime = time.time() - start
elapsedTime = calcElapsedTime(currentTime)
print(f"Elapsed time: {elapsedTime}")
trackerOut = open(os.path.join(workingDir, outputPrefix + "_sppIDerRun.info"), 'a')
trackerOut.write(f"\nMQ scores plotted\nElapsed time: {elapsedTime}")
trackerOut.close()

########################## bedgraph Coverage ###########################
sortOut = bamSortOut
if args.bed:
    bedOutD = outputPrefix + "-d.bedgraph"
    bedFileD = open(os.path.join(workingDir, bedOutD), 'w')
    subprocess.call(["genomeCoverageBed", "-d", "-ibam", sortOut], stdout=bedFileD, cwd=workingDir)
    bedFileD.close()
else:
    bedOut = outputPrefix + ".bedgraph"
    bedFile = open(os.path.join(workingDir, bedOut), 'w')
    subprocess.call(["genomeCoverageBed", "-bga", "-ibam", sortOut], stdout=bedFile, cwd=workingDir)
    bedFile.close()
print("bedgraph complete")
currentTime = time.time() - start
elapsedTime = calcElapsedTime(currentTime)
print(f"Elapsed time: {elapsedTime}")
trackerOut = open(os.path.join(workingDir, outputPrefix + "_sppIDerRun.info"), 'a')
trackerOut.write(f"\nbedgraph complete\nElapsed time: {elapsedTime}")
trackerOut.close()

########################## average Bed ###########################
if args.bed:
    subprocess.call(["Rscript", os.path.join(scriptDir, "meanDepth_sppIDer-d.R"), outputPrefix], cwd=workingDir)
else:
    subprocess.call(["Rscript", os.path.join(scriptDir, "meanDepth_sppIDer-bga.R"), outputPrefix], cwd=workingDir)
print("Found mean depth")
currentTime = time.time() - start
elapsedTime = calcElapsedTime(currentTime)
print(f"Elapsed time: {elapsedTime}")
trackerOut = open(os.path.join(workingDir, outputPrefix + "_sppIDerRun.info"), 'a')
trackerOut.write(f"\nFound mean depth\nElapsed time: {elapsedTime}")
trackerOut.close()

########################## make plot ###########################
subprocess.call(["Rscript", os.path.join(scriptDir, "sppIDer_depthPlot_forSpc.R"), outputPrefix], cwd=workingDir)
subprocess.call(["Rscript", os.path.join(scriptDir, "sppIDer_depthPlot.R"), outputPrefix], cwd=workingDir)
print("Plot complete")
currentTime = time.time() - start
elapsedTime = calcElapsedTime(currentTime)
print(f"Elapsed time: {elapsedTime}")
trackerOut = open(os.path.join(workingDir, outputPrefix + "_sppIDerRun.info"), 'a')
trackerOut.write(f"\nPlot complete\nElapsed time: {elapsedTime}\n")
trackerOut.close()

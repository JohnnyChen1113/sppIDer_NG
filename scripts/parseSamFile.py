__author__ = 'Quinn'
__modified_by__ = 'Junhao Chen'

import sys, re, time, os
from pathlib import Path

################################################################
# This script will take the sam formatted output of bwa and parse for mapping quality and which genomes reads map to.
#
# Input: sam file of reads mapped to a combination reference genome
################################################################

# 获取脚本所在目录和当前工作目录
scriptDir = Path(__file__).resolve().parent
workingDir = Path.cwd()

inputName = sys.argv[1]
samName = inputName + ".sam"
outputName = inputName + "_MQ.txt"
outputLenName = inputName + "_chrLens.txt"
start = time.time()

speciesDict = {}
MQscoreDict = {}
speciesDict["*"] = {}
speciesDict["*"][0] = 0
speciesList = ['*']

with open(workingDir / outputLenName, 'w', encoding='utf-8') as outputLen, open(workingDir / samName, 'r', encoding='utf-8') as sam:
    samLines = sam.read().splitlines()
    for line in samLines:
        if re.match('^(@SQ)', line):
            headerInfo = line.split('\t')
            chrInfo = headerInfo[1].split(":")[1]
            chrName = chrInfo.split("-")
            speciesName = chrName[0]
            chrNum = int(chrName[1])
            chrLen = headerInfo[2].split(":")[1]
            outputLen.write(chrInfo + "\t" + chrLen + "\n")
            if chrNum == 1:
                speciesList.append(speciesName)
                speciesDict[speciesName] = {}
                for i in range(0, 61):
                    speciesDict[speciesName][i] = {'count': 0, 'names': []}
        elif re.match('^(?!@)', line):
            lineSplit = line.split('\t')
            sequenceName = lineSplit[0]
            chr = lineSplit[2].split("-")
            species = chr[0]
            MQscore = int(lineSplit[4])
            speciesDict[species][MQscore]['count'] += 1
            speciesDict[species][MQscore]['names'].append(sequenceName)

with open(workingDir / outputName, 'w', encoding='utf-8') as output:
    output.write("Species\tMQscore\tcount\tSequenceNames\n")
    for species in speciesList:
        for score in speciesDict[species].keys():
            count = speciesDict[species][score]['count']
            names = ",".join(speciesDict[species][score]['names'])
            output.write(f"{species}\t{score}\t{count}\t{names}\n")

currentTime = time.time() - start
print(str(currentTime) + " secs\n")


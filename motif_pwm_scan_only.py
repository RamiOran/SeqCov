# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
from __future__ import division
"""
Perform operations to ONLY scan PWM on input file   and find how much coverage by all the motifs in the PWM file
"""
import Fimo
import general_utils
import re


def getPwmOrder(pwmFileName):
	"""
	Read the PWm file and get the motifs order as in the file
	"""
	motifList = []
	with open(pwmFileName, 'rb') as handler:
		for line in handler:
			if re.search(r'MOTIF', line):
				split = line.split()
				motifId = split[1]
				motifList.append(motifId)
	
	print 'motifList:', motifList	
	return motifList


def callMotifPwmScanOnly(args, confDict, fileList):
	"""
	
	"""
	#get the PWM file to scan for 
	pwmFileName = confDict['input']['pwm_file']
	#find the total number of seqs in the testing file
	foreNumSeqs = general_utils.findNumSeqs(confDict['input']['fore_testing_file']) 
	if confDict['input']['back_testing_file'] != 'none':
		backNumSeqs = general_utils.findNumSeqs(confDict['input']['back_testing_file'])
	
	#read the PWM file and get the order of the motifs so to be used when reporting the results
	orderMotifList = getPwmOrder(pwmFileName)
	
	#run the FIMO tool for froeground file
	fimoParaDict = {}
	fimoParaDict['thresh'] = confDict['fimo']['pvalue']
	#fimo output directory name
	fimoForeOutDir = args.jid + '_fore_Fimo'
	Fimo.callFimo(confDict['input']['fore_testing_file'] , pwmFileName, fimoParaDict, fimoForeOutDir)
	fileList.append(fimoForeOutDir)
	#parse the fimo results folder
	strand = confDict['fimo']['strand']
	#fimoDict between motif ID and seq hits in the testing file
	foreFimoDict = Fimo.parseFimo(fimoForeOutDir+'/fimo.txt', strand)
	#write the motifs in the format for sequence coverge algorithms 
	foreMotifSeqFileName = args.jid + '_fore_motif_hits'
	general_utils.writeMotifSeqFile(foreFimoDict, foreMotifSeqFileName)
	fileList.append(foreMotifSeqFileName)
	
	
	#run the FIMO tool for the background file if it exists
	if confDict['input']['back_testing_file'] != 'none':
		#Background fimo output dir name
		fimoBackOutDir = args.jid + '_back_Fimo'
		Fimo.callFimo(confDict['input']['back_testing_file'], pwmFileName, fimoParaDict, fimoBackOutDir)
		fileList.append(fimoBackOutDir)
		#parse the fimo results folder
		strand = confDict['fimo']['strand']
		#fimoDict between motif ID and seq hits in background file
		backFimoDict = Fimo.parseFimo(fimoBackOutDir + '/fimo.txt', strand)
		#write the motifs in the format for sequence covearge algorithms
		backMotifSeqFileName = args.jid + '_back_motif_hits'
		general_utils.writeMotifSeqFile(backFimoDict, backMotifSeqFileName)
		fileList.append(backMotifSeqFileName) 
	
	
	outFileName = args.jid + '_scan_results.csv'
	fileList.append(outFileName)
	outFile = open(outFileName, 'wb')
	if confDict['input']['back_testing_file'] != 'none':
		outFile.write('#Name,foreground_motif_seqCount,foreground_motif_seqCoverage(%),foreground_seqs_added,cumulative_foreground_seq_cov(%)' 
		   + ',background_motif_seqCount,background_motif_seqCoverage(%),background_seqs_added,cumulative_background_seq_cov(%)' +'\n')
	else:
		outFile.write('#Name,foreground_motif_seqCount,foreground_motif_seqCoverage(%),foreground_seqs_added,cumulative_foreground_seq_cov(%)'+'\n')
	
	#go thru the fimo directory and write the coverage information
	counter = 1
	allSeqSet = set()
	backAllSeqSet = set()
	backCumCov = 0
	#for motifId in foreFimoDict:
	for motifId in orderMotifList:
		lineList = []
		print 'm:', motifId
		if motifId not in foreFimoDict:
			continue
		#foreground file
		lineList.append(motifId)
		seqList = foreFimoDict[motifId]
		seqCount = len(seqList)
		lineList.append(str(seqCount))
		cov = 100*(seqCount/foreNumSeqs)
		lineList.append(str(cov))
		seqSet = set(seqList)
		seqDiff =  seqSet - allSeqSet
		newSeqsAdded = len(seqDiff)
		allSeqSet = allSeqSet | seqSet#diff in sets
		cumCov = 100*(len(allSeqSet)/foreNumSeqs)
		lineList.append(str(newSeqsAdded))
		lineList.append(str(cumCov))
		#background information
		if confDict['input']['back_testing_file'] != 'none':
			if motifId in backFimoDict:
				backSeqList = backFimoDict[motifId]
				backSeqCount = len(backSeqList)
				lineList.append(str(backSeqCount))
				backCov = 100*(backSeqCount/backNumSeqs)
				lineList.append(str(backCov))
				backSeqSet = set(backSeqList)
				backSeqDiff = backSeqSet - backAllSeqSet
				backNewSeqsAdded = len(backSeqDiff)
				backAllSeqSet = backAllSeqSet | backSeqSet
				backCumCov = 100*(len(backAllSeqSet)/backNumSeqs)
				lineList.append(str(backNewSeqsAdded))
				lineList.append(str(backCumCov))
			else:
				backSeqCount = 0
				lineList.append(str(backSeqCount))
				backCov = 0
				lineList.append(str(backCov))
				backNewSeqsAdded = 0
				backCumCov += 0
				lineList.append(str(backNewSeqsAdded))
				lineList.append(str(backCumCov))
				
		line = ','.join(lineList)
		
		outFile.write(line + '\n')
		counter += 1

	outFile.close()


def main():
	pass


if( __name__ == "__main__" ):
    main(sys.argv)
    pass
else:
    pass

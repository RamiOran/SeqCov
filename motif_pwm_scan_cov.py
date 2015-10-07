# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
from __future__ import division
"""
Perform operations to scan PWM input file and find coverage and number of motifs required (The motif selection problem)
"""
#python imports
import sys
import os 
import argparse
import shutil

#add the util folder path to use util files in it
inPath = os.path.realpath(__file__)
split = inPath.split('/')
inPath = '/'.join(split[:len(split)-1])
sys.path.append(inPath + '/utils')
import conf
import general_utils
import Fimo
import Tomtom
#add the alg folder path to call the different algorithms
sys.path.append('./algo')
import greedy
import ILP
import branch_cut
import required_cover


def outputResults(outFileName, depthDict, backFimoDict, filterMinNumSeq, idMotifDict, foreNumSeqs, backNumSeqs):
	"""Output the results from the depth dictionary
	Args:
		outFileName: name of file to output results to
	"""
	# a list that stores the fianl list of IDs that passed the threshold filter
	finalSelectMotifList = []
	depthOut = open(outFileName, 'wb')
	for depthVal in sorted(depthDict.iterkeys(), key=int):
		print 'depth is:', depthVal
		depthOut.write('\nDepth:' + str(depthVal) + '\n')
		depthOut.write('#Name,foreground_motif_seqCount,foreground_motif_seqCoverage(%),foreground_seqs_added,cumulative_foreground_seq_cov(%)' 
		   + ',background_motif_seqCount,background_motif_seqCoverage(%),background_seqs_added,cumulative_background_seq_cov(%)' +'\n')
		cumCov = 0
		backCumCov = 0
		backSeqsAdded = []
		motifIdList = depthDict[depthVal][0]
		motifSeqSet = depthDict[depthVal][1]
		newSeqsAddedList = depthDict[depthVal][2]
		#add new information to the motif objects as well
		for i in range(len(motifIdList)):
			lineList = []
			motifId = motifIdList[i]
			motifSeqset = motifSeqSet[i]
			seqsAdded = newSeqsAddedList[i]
			#apply the filter threshold
			if seqsAdded < filterMinNumSeq:
				continue
			if motifId not in idMotifDict:
				motifName = motifId
			else:
				motifName = idMotifDict[motifId]
			finalSelectMotifList.append(motifName)
			
			lineList.append(motifName)
			lineList.append(str(len(motifSeqset)))
			cov = 100 * (len(motifSeqset)/foreNumSeqs)
			lineList.append(str(cov))
			lineList.append(str(seqsAdded))
			cumCov += (seqsAdded/foreNumSeqs)*100
			lineList.append(str(cumCov))
			#add the background information for this motif
			#in the backFimoDict the IDs are string type
			
			if motifName in backFimoDict:
				backSeqList = backFimoDict[motifName]
				backSeqCount = len(backSeqList)
				backSeqCov = (backSeqCount/backNumSeqs)*100
				backSeqsAddedCount = 0
				for backSeq in backSeqList:
					if backSeq not in backSeqsAdded:
						backSeqsAddedCount += 1
						backSeqsAdded.append(backSeq)
				backCumCov += (backSeqsAddedCount/backNumSeqs)*100
			else:
				backSeqCount = 0
				backSeqCov = 0
				backSeqsAddedCount=0
				backCumCov += 0
			lineList.append(str(backSeqCount))
			lineList.append(str(backSeqCov))
			lineList.append(str(backSeqsAddedCount))
			lineList.append(str(backCumCov))
			line = ','.join(lineList)
			depthOut.write(line + '\n')
	depthOut.close()
	return finalSelectMotifList


def callMotifPwmScanCov(args, confDict, fileList):
	"""
	The main function to perfrom motif PWM scanning and coverage on a testing file
	Args:
		args: the object that holds the input args
		confDict: the dict of the input configuration file
		fileList: list of files to be added to the results directory
	Returns:
	"""
	#get the PWM file to scan for 
	pwmFileName = confDict['input']['pwm_file']
	#read the PWM file and make a motif dict between motif IDs (names) and the motif objects which have PWM lines. Also return two dicts between the motif anmes and IDs and vice versa
	motifDict = general_utils.makeMotifDictFromPwmFile(pwmFileName)
	print 'len motifDict:', len(motifDict)
	#find the total number of seqs in the testing file
	foreNumSeqs = general_utils.findNumSeqs(confDict['input']['fore_testing_file'])
	#make a list of all seq names in the foreground data set
	foreSeqList = general_utils.findSeqList(confDict['input']['fore_testing_file'])
	#make a list of all seq names in the background data set
	backSeqList = general_utils.findSeqList(confDict['input']['back_testing_file'])
	
	#run the FIMO tool on the foreground data set
	fimoParaDict = {}
	fimoParaDict['thresh'] = confDict['fimo']['pvalue']
	#fimo output directory name
	fimoOutDir = args.jid + '_foreground_Fimo'
	Fimo.callFimo(confDict['input']['fore_testing_file'] , pwmFileName, fimoParaDict, fimoOutDir)
	fileList.append(fimoOutDir)
	#parse the fimo results folder
	strand = confDict['fimo']['strand']
	#fimoDict between motif ID and seq hits in the foeground testing file
	fimoDict = Fimo.parseFimo(fimoOutDir+'/fimo.txt', strand)
	#write the motifs in the format for sequence coverge algorithms 
	motifSeqFileName = args.jid + '_fore_motif_hits_in_seqs'
	general_utils.writeMotifSeqFile(fimoDict, motifSeqFileName)
	fileList.append(motifSeqFileName)
	
	
	#check if motif refinement requested to remove redundant motifs when applying the seq_cov algthms
	#the tomtomDict will be used with the seq coverage algorithm. The dict is between all the motifs to check their similarity 
	tomtomDict = {}
	if confDict['motif_refine']['app'] == 'tomtom_refine':
		print 'Peform motif refine using tomtom'
		#dict that stores parameters for the tomtom tool
		tomParaDict = {}
		tomParaDict['evalue'] = confDict['tomtom_refine']['evalue']
		#tomtom output directory name
		tomtomOutDir = args.jid + '_Tomtom_refine'
		#we are comparing the motifs against each other, so use the same PWM file
		Tomtom.callTomtom(confDict['input']['pwm_file'], confDict['input']['pwm_file'], tomParaDict, tomtomOutDir)
		fileList.append(tomtomOutDir)
		#parse the tomtom output file and return a dict between motif Name and list[] of motifs that match it
		tomtomDict = Tomtom.parse(tomtomOutDir+'/tomtom.txt')
	
	
	#process the seq cov file
	motifIdDict, idMotifDict, seqIdDict, idSeqDict, Uset, Sdict = general_utils.processMotifHitFile_1(motifSeqFileName)
	#call the seq cov algthm
	minSetIdList = [] 
	minSeqSet = set()
	
	
	#set the filtering values for later
	filterThresh = 0
	filterMinNumSeq = 0	
	
	#this dict is between a depth value and a tuple of motif Ids, seq sets, and newly added seqs by each motif
	depthDict = {}
	
	
	#just a tmp check will be removed later
	makeLogo = 'false'
	
	
	#####
	#ILP 
	####
	if confDict['sequence.coverage']['app'] == 'ILP' and confDict['motif_refine']['app'] == 'none':
		print 'ILP algorithm'
		depth = int(confDict['ILP']['depth'])
		#this dict is between a depth value and a list of motif IDs
		depthDict = ILP.callILPDepth(args.jid, fimoDict, tomtomDict, depth, fileList)
		print 'len:', len(depthDict)
		#if there is a background file use it for comparison purposes
		if confDict['input']['back_testing_file'] != 'none':
			#run the FIMO tool for the background sequences
			fimoParaDict = {}
			fimoParaDict['thresh'] = confDict['fimo']['pvalue']
			#fimo output directory name
			fimoOutDir = args.jid + '_background_Fimo'
			Fimo.callFimo(confDict['input']['back_testing_file'] , pwmFileName, fimoParaDict, fimoOutDir)
			fileList.append(fimoOutDir)
			#find number of sequences in the background file
			backNumSeqs = general_utils.findNumSeqs(confDict['input']['back_testing_file'])
			#parse the fimo results folder
			strand = confDict['fimo']['strand']
			#fimoDict between motif names and seq hits in the background  file
			backFimoDict = Fimo.parseFimo(fimoOutDir+'/fimo.txt', strand)
			#write the motifs in the format for sequence coverge algorithms 
			motifSeqFileName = args.jid + '_back_motif_hits_in_seqs'
			general_utils.writeMotifSeqFile(backFimoDict, motifSeqFileName)
			fileList.append(motifSeqFileName)
		
		#no filtering
		filterMinNumSeq = 0
		#output results 
		for depthVal in sorted(depthDict.iterkeys(), key=int):
			depthOutFileName = args.jid + '_depth_results_'+ str(depthVal) + '.csv'
			finalSelectMotifList = outputResults(depthOutFileName, depthDict, backFimoDict, filterMinNumSeq, idMotifDict,foreNumSeqs, backNumSeqs)
			print 'final:', finalSelectMotifList
			fileList.append(depthOutFileName)
		
	####	
	#branch and cut	
	####
	if confDict['sequence.coverage']['app'] == 'branch_cut' and confDict['motif_refine']['app'] == 'none':
		print 'Branch and cut algorithm'
		depth = int(confDict['branch_cut']['depth'])
		#this dict is between a depth value and a list of motif IDs
		depthDict = branch_cut.callBranchDepth(args.jid, fimoDict, tomtomDict, depth, fileList)
		print 'len:', len(depthDict)
		
		#if there is a background file use it for comparison purposes
		if confDict['input']['back_testing_file'] != 'none':
			#run the FIMO tool for the background sequences
			fimoParaDict = {}
			fimoParaDict['thresh'] = confDict['fimo']['pvalue']
			#fimo output directory name
			fimoOutDir = args.jid + '_background_Fimo'
			Fimo.callFimo(confDict['input']['back_testing_file'] , pwmFileName, fimoParaDict, fimoOutDir)
			fileList.append(fimoOutDir)
			#find number of sequences in the background file
			backNumSeqs = general_utils.findNumSeqs(confDict['input']['back_testing_file'])
			#parse the fimo results folder
			strand = confDict['fimo']['strand']
			#fimoDict between motif names and seq hits in the background  file
			backFimoDict = Fimo.parseFimo(fimoOutDir+'/fimo.txt', strand)
			#write the motifs in the format for sequence coverge algorithms 
			motifSeqFileName = args.jid + '_back_motif_hits_in_seqs'
			general_utils.writeMotifSeqFile(backFimoDict, motifSeqFileName)
			fileList.append(motifSeqFileName)
		
		#no filtering
		filterMinNumSeq = 0
		#output results 
		for depthVal in sorted(depthDict.iterkeys(), key=int):
			depthOutFileName = args.jid + '_depth_results_'+ str(depthVal) + '.csv'
			finalSelectMotifList = outputResults(depthOutFileName, depthDict, backFimoDict, filterMinNumSeq, idMotifDict,foreNumSeqs, backNumSeqs)
			print 'final:', finalSelectMotifList
			fileList.append(depthOutFileName)
	
	####	
	#required cover	
	####
	if confDict['sequence.coverage']['app'] == 'required_cover' and confDict['motif_refine']['app'] == 'none':
		print 'Required cover algorithm'
		depth = int(confDict['required_cover']['depth'])
		#this dict is between a depth value and a list of motif IDs
		depthDict = required_cover.callRequiredDepth(args.jid, fimoDict, tomtomDict, depth, fileList)
		print 'len:', len(depthDict)
		#if there is a background file use it for comparison purposes
		if confDict['input']['back_testing_file'] != 'none':
			#run the FIMO tool for the background sequences
			fimoParaDict = {}
			fimoParaDict['thresh'] = confDict['fimo']['pvalue']
			#fimo output directory name
			fimoOutDir = args.jid + '_background_Fimo'
			Fimo.callFimo(confDict['input']['back_testing_file'] , pwmFileName, fimoParaDict, fimoOutDir)
			fileList.append(fimoOutDir)
			#find number of sequences in the background file
			backNumSeqs = general_utils.findNumSeqs(confDict['input']['back_testing_file'])
			#parse the fimo results folder
			strand = confDict['fimo']['strand']
			#fimoDict between motif names and seq hits in the background  file
			backFimoDict = Fimo.parseFimo(fimoOutDir+'/fimo.txt', strand)
			#write the motifs in the format for sequence coverge algorithms 
			motifSeqFileName = args.jid + '_back_motif_hits_in_seqs'
			general_utils.writeMotifSeqFile(backFimoDict, motifSeqFileName)
			fileList.append(motifSeqFileName)
		
		#no filtering
		filterMinNumSeq = 0
		#output results 
		for depthVal in sorted(depthDict.iterkeys(), key=int):
			depthOutFileName = args.jid + '_depth_results_'+ str(depthVal) + '.csv'
			finalSelectMotifList = outputResults(depthOutFileName, depthDict, backFimoDict, filterMinNumSeq, idMotifDict,foreNumSeqs, backNumSeqs)
			print 'final:', finalSelectMotifList
			fileList.append(depthOutFileName)
		
	
	####
	#greedy filtered (greedyFilt) and refinement	
	####
	if confDict['sequence.coverage']['app'] == 'greedyFilt' and confDict['motif_refine']['app'] == 'tomtom_refine':
		depth = int(confDict['greedyFilt']['depth'])
		if confDict['greedyFilt']['filter'] == 'true':
			#get the threshold value from the configure file
			filterThresh = float(confDict['greedyFilt']['filter_threshold'])
			filterMinNumSeq = filterThresh * foreNumSeqs
		print 'Perform filtered greedy sequence coverage with refinement and depth:', depth
		#this dict is between a depth value and a tuple of motif Ids, seq sets, and newly added seqs by each motif
		#depthDict = greedy.callGreedyDepthRefine(Uset, Sdict, depth, tomtomDict, motifIdDict, idMotifDict, filterMinNumSeq)
		depthDict = greedy.callGreedyFiltRefine(Uset, Sdict, depth, tomtomDict, motifIdDict, idMotifDict, filterMinNumSeq)
		#if there is a background file use it for comaprison purposes
		if confDict['input']['back_testing_file'] != 'none':
			#run the FIMO tool for the background sequences
			fimoParaDict = {}
			fimoParaDict['thresh'] = confDict['fimo']['pvalue']
			#fimo output directory name
			fimoOutDir = args.jid + '_background_Fimo'
			Fimo.callFimo(confDict['input']['back_testing_file'] , pwmFileName, fimoParaDict, fimoOutDir)
			fileList.append(fimoOutDir)
			#find number of sequences in the background file
			backNumSeqs = general_utils.findNumSeqs(confDict['input']['back_testing_file'])
			#parse the fimo results folder
			strand = confDict['fimo']['strand']
			#fimoDict between motif names and seq hits in the background  file
			backFimoDict = Fimo.parseFimo(fimoOutDir+'/fimo.txt', strand)
			#write the motifs in the format for sequence coverge algorithms 
			motifSeqFileName = args.jid + '_back_motif_hits_in_seqs'
			general_utils.writeMotifSeqFile(backFimoDict, motifSeqFileName)
			fileList.append(motifSeqFileName)
			
		
			
		#output results 
		depthOutFileName = args.jid + '_depth_results.csv'
		finalSelectMotifList = outputResults(depthOutFileName, depthDict, backFimoDict, filterMinNumSeq, idMotifDict, foreNumSeqs, backNumSeqs)
		fileList.append(depthOutFileName)
	
	
	####
	#check if want to produce motif logos
	####
	if confDict['motif.logo']['opt'] == 'true':
		if makeLogo == 'false':	
			#get all the motif IDs from all the depths
			mIdList = []
			for depthVal in sorted(depthDict.iterkeys(), key=int):
				mIdList.extend(depthDict[depthVal][0]) 
			seqCovPwmFileName = args.jid + '_seqCovMotifs.pwm'
			logoDirName = args.jid + '_seqCovLogo'
			fileList.append(seqCovPwmFileName)
			fileList.append(logoDirName)
			general_utils.makeMotifLogo(mIdList, motifDict, logoDirName, seqCovPwmFileName, idMotifDict, 'cov', finalSelectMotifList)	
	
			
def main():
	pass


if( __name__ == "__main__" ):
    main(sys.argv)
else:
    pass

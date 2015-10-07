# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
""" Greedy algortihm for motif selection"""

from __future__ import division
import sys
from copy import deepcopy
##



def findMax(Sdict, R, foundIdList):
	"""Find the motif(subset) that covers the most uncovered sequences.
	
	Args:
		Sdict: dict between motif Ids and sequences the motif occurs in
		R: the list of sequences not covered yet
		foundIdList: list of motifs simialr to previous motifs reported
		
	Returns:
		subset with which covers the most and motif ID of this subset
	"""
	maxCov = 0 
	newSeqsAdded = 0
	maxCovMotif = -1
	for motifId, sSet in Sdict.items():
		#skip if motif already reported or is simialr to previous motifs
		if motifId in foundIdList:
			continue
		if len(sSet.intersection(R)) == 0:#no new 
			continue
		else:
			#print '\tId:', motifId, sSet
			#tmpSet = sSet.intersection(R)
			cov = len(sSet.intersection(R))
			if cov > maxCov:
				maxCov = cov
				maxCovMotif = motifId
				newSeqsAdded = len(sSet.intersection(R))  
	
	if maxCovMotif == -1:
		return list(), -1, list()
	#print 'miElem:', minElement,'Sdict[min]', Sdict[minElement]
	#Sdict[minElement] is the list of seqs covered by motif maxElement 
	return Sdict[maxCovMotif], maxCovMotif, newSeqsAdded
	
def callGreedyFiltRefine(Uset, Sdict, depth, tomtomDict, motifIdDict, idMotifDict, filterMinNumSeq):
	"""Call the greedy method where the filter step is integrated into the algortihm itself """
	print 'Filtered greedy sequence coverage'
	#see how many motifs (subsets) we have
	numMotifs = len(Sdict)
	#make a list of weigts equal to zero
	wghtList  = [1] * numMotifs
	print 'filterMinNumSeq:', filterMinNumSeq
	R = Uset
	C = []
	j = 0
	motifIdList = []
	newSeqsAddedList = []
	foundIdList = []
	while len(R) != 0 and j < len(Sdict):
		S_i, motifId, newSeqsAdded = findMax(Sdict, R, foundIdList)
		if motifId != -1:
			if newSeqsAdded < filterMinNumSeq:
				break
			#remove similar motifs to this one 
			for val in tomtomDict:
				if isinstance(val,str):
					tomIdType = 'str'
				if isinstance(val,int):
					tomIdType = 'integer'
				break
				
			motifName = idMotifDict[motifId]
			for matchName in tomtomDict[motifName]:
				if matchName not in motifIdDict:
					continue
				matchId = motifIdDict[matchName]
				foundIdList.append(matchId)
			
			
			C.append(S_i)
			motifIdList.append(motifId)
			R = R.difference(S_i)
			newSeqsAddedList.append(newSeqsAdded)
			j += 1
		else:
			break
	
	#print 'Cover:', C
	print 'motif ID list:', motifIdList
	
	#print "Total Cost: ", sum(costs), costs
	#return motifIdList, C
	#newSeqsAddedList =[200] * len(motifIdList)
	depthDict = {}
	print motifIdList
	for val in motifIdList:
		print val,idMotifDict[val]
	depthDict[1] = (motifIdList, C, newSeqsAddedList)
	return depthDict



def main(args):
    pass

if( __name__ == "__main__" ):
    main(sys.argv)
    pass
else:
    pass


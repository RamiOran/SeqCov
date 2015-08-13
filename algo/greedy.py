from __future__ import division
"""
greedy.py -- greedy algorithm
greedy.py -- greedy algorithm

Changelog:
    greedy.py ; started at 2014-09-24 16:44:57 EDT by Liang Chen
"""
import sys
from copy import deepcopy
##

def findMin(Sdict, R, w):
	"""Find the min cost subset. the subset that covers the most seqiences is the best.
	
	Args:
		Sdict: dict between motif Ids and sequences the motif occurs in
		R: the list of sequences not covered yet
		w: list of weights assigned to each subset/motif. For greedy here all set to 1 (unicost set cover)
	
	Returns:
		subset with lowest cost and cost of it and motif ID of this subset
	
	"""
	minCost = 99999 
	minElement = -1
	for motifId, sSet in Sdict.items():
		if len(sSet.intersection(R)) == 0:#no new 
			continue
		else:
			#print '\tId:', motifId, sSet
			tmpSet = sSet.intersection(R)
			#print '\ttmpSet:', tmpSet
			cost = w[motifId-1]/(len(sSet.intersection(R)))
			if cost < minCost:
				minCost = cost
				minElement = motifId 
	
	#print 'miElem:', minElement,'Sdict[min]', Sdict[minElement]
	return Sdict[minElement], w[minElement-1], minElement
	


def callGreedy(Uset, Sdict):
	"""use the greedy algthm to find the smallest set cover
	
	Follow the pseudo code by the introduction to algorithms book and the wiki page. Very simple
	Just pick the set which covers the most at each iteration.
	Used an online implementation for now. This is greedy with no weights (unicost set cover) so make the weights
	vector all equal to one
	
	Args:
		Uset: set of seq Ids
		Sdict; dict between motif Ids and the set of sequences each motif
		
	Returns:
		motifIdList: list of motif Ids with smallest set coverage
		C: the seq sets found, 
	"""
	print "greedy.py::call_greedy()"
	#see how many motifs (subsets) we have
	numMotifs = len(Sdict)
	#make a list of weigts equal to zero
	wghtList  = [1] * numMotifs
    
   
	R = Uset
	C = []
	costs = []
	counter = 0
	motifIdList = []
	while len(R) != 0:
		S_i, cost, motifId = findMin(Sdict, R, wghtList)
		C.append(S_i)
		motifIdList.append(motifId)
		R = R.difference(S_i)
		costs.append(cost)
		counter += 1
	
	#print 'Cover:', C
	print 'motif ID list:', motifIdList
	#print "Total Cost: ", sum(costs), costs
	return motifIdList, C
    

def findMinDepth(Sdict, R, w, foundIdList):
	"""Find the min cost subset. the subset that covers the most seqiences is the best.
	This has a new list attached to it which is the list of IDs already found
	Args:
		Sdict: dict between motif Ids and set of sequences the motif occurs in
		R: the list of sequences not covered yet
		w: list of weights assigned to each subset/motif. For greedy here all set to 1 (unicost set cover)
	
	Returns:
		subset with lowest cost and cost of it and motif ID of this subset
	
	"""
	minCost = 99999 
	minElement = -1
	newSeqsAdded = 0
	for motifId, sSet in Sdict.items():
		#skip if motif already reported
		if motifId in foundIdList:
			continue
		if len(sSet.intersection(R)) == 0:#no new 
			continue
		else:
			#print '\tId:', motifId, sSet
			#tmpSet = sSet.intersection(R)
			cost = w[motifId-1]/(len(sSet.intersection(R)))
			if cost < minCost:
				minCost = cost
				minElement = motifId
				newSeqsAdded = len(sSet.intersection(R)) 
	
	if minElement == -1:
		return list(), list(), -1, list()
	#print 'miElem:', minElement,'Sdict[min]', Sdict[minElement]
	return Sdict[minElement], w[minElement-1], minElement, newSeqsAdded


def callGreedyDepth(Uset, Sdict, depth):
	"""use the greedy algthm to find the smallest set cover with depth. For each depth discovered remove the motif IDs found
	for that dpeth then do the search again
	
	Args:
		Uset: set of seq Ids
		Sdict; dict between motif Ids and the set of sequences each motif
		
	Returns:
		motifIdList: list of motif Ids with smallest set coverage
		C: the seq sets found, 
	"""
	print 'Greedy sequence coverage with depth:', depth
	#see how many motifs (subsets) we have
	numMotifs = len(Sdict)
	#make a list of weigts equal to zero
	wghtList  = [1] * numMotifs
	#make a list that stores motif IDs found so far
	foundIdList = []
	#make a dpth dictionary between the depth and a tuple of (motifIdList, list of sets of seqs per motif)
	depthDict = {}
    #loop thru the depths 
	for i in range(1,depth+1):
		print 'depth:', i
		R = Uset
		C = []
		costs = []
		counter = 0
		motifIdList = []
		newSeqsAddedList = []
		while len(R) != 0:
			S_i, cost, motifId, newSeqsAdded = findMinDepth(Sdict, R, wghtList,foundIdList)
			C.append(S_i)
			motifIdList.append(motifId)
			R = R.difference(S_i)
			newSeqsAddedList.append(newSeqsAdded)
			costs.append(cost)
			counter += 1
		
		#print 'Cover:', C
		print 'motif ID list:', motifIdList
		#update the found list
		for newId in motifIdList:
			foundIdList.append(newId)
		#make a tuple and store into depth dict
		depthDict[i] = (motifIdList, C, newSeqsAddedList)
		#print "Total Cost: ", sum(costs), costs
		#return motifIdList, C
	
	#return the depth dictionary
	return depthDict




def findMinDepthRefine(Sdict, R, w, foundIdList):
	"""Find the min cost subset. the subset that covers the most seqiences is the best.
	This has a new list attached to it which is the list of IDs already found
	Args:
		Sdict: dict between motif Ids and sequences the motif occurs in
		R: the list of sequences not covered yet
		w: list of weights assigned to each subset/motif. For greedy here all set to 1 (unicost set cover)
	
	Returns:
		subset with lowest cost and cost of it and motif ID of this subset
	
	"""
	minCost = 99999 
	minElement = -1
	newSeqsAdded = 0
	for motifId, sSet in Sdict.items():
		#skip if motif already reported or is simialr to previous motifs
		if motifId in foundIdList:
			continue
		if len(sSet.intersection(R)) == 0:#no new 
			continue
		else:
			#print '\tId:', motifId, sSet
			#tmpSet = sSet.intersection(R)
			cost = w[motifId-1]/(len(sSet.intersection(R)))
			if cost < minCost:
				minCost = cost
				minElement = motifId
				newSeqsAdded = len(sSet.intersection(R))  
	
	if minElement == -1:
		return list(), list(), -1, list()
	#print 'miElem:', minElement,'Sdict[min]', Sdict[minElement]
	return Sdict[minElement], w[minElement-1], minElement, newSeqsAdded


def callGreedyDepthRefine(Uset, Sdict, depth, tomtomDict, motifIdDict, idMotifDict, filterMinNumSeq):
	"""use the greedy algthm to find the smallest set cover with depth. For each depth discovered remove the motif IDs found
	for that dpeth then do the search again
	Here with motif refinement applied such that simialr motifs are removed 
	Args:
		Uset: set of seq Ids
		Sdict: dict between motif Ids and the set of sequences each motif
		depth: coverage depth
		tomtomDict: dict between motif Ids and a list of motif Ids that are similar to it
		motifIdDict: dict between motif name and ID (as found by processMotifHitFile)
		idMotifDict: dict between motif id and motif name (as found by processMotifHitFile)
	Returns:
		depthDict: dictionary between the depth and a tuple of (motifIdList, list of sets of seqs per motif)
	"""
	print 'Greedy sequence coverage with depth:', depth
	#see how many motifs (subsets) we have
	numMotifs = len(Sdict)
	#make a list of weigts equal to 1
	wghtList  = [1] * numMotifs
	#make a depth dictionary between the depth and a tuple of (motifIdList, list of sets of seqs per motif)
	depthDict = {}
	#ID of motifs found per depth
	motifIdList = []
    #loop thru the depths 
	for i in range(1,depth+1):
		print 'depth:', i
		#make a list that stores motif IDs found so far for depth coverage
		foundIdList = []
		#update the found list
		for newId in motifIdList:
			foundIdList.append(newId)
		print 'foundIdList:', foundIdList
		motifIdList = []
		R = Uset
		C = []
		costs = []
		newSeqsAddedList = []
		while len(R) != 0:
			#S_i is a set
			S_i, cost, motifId, newSeqsAdded = findMinDepthRefine(Sdict, R, wghtList, foundIdList)
			#check if nothing found so go to next depth
			if motifId != -1:
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
				
				#check how many new seqs added
				print 'newAdded:', newSeqsAdded,'filter:', filterMinNumSeq
				if newSeqsAdded > filterMinNumSeq:
					C.append(S_i)
					newSeqsAddedList.append(newSeqsAdded)
					R = R.difference(S_i)
					costs.append(cost)
					motifIdList.append(motifId)
				else:
					break
					
			else:
				break
		
		#print 'Cover:', C
		#make a tuple and store into depth dict
		print 'after list:', motifIdList
		depthDict[i] = (motifIdList, C, newSeqsAddedList)
		#print "Total Cost: ", sum(costs), costs
		#return motifIdList, C
	
	#return the depth dictionary
	return depthDict


def main(args):
    pass

if( __name__ == "__main__" ):
    main(sys.argv)
    pass
else:
    pass


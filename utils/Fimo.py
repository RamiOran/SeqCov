# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
from __future__ import division
"""
Used to call the FIMO tool; part of the MEME suite
http://meme-suite.org/
"""
import os
import re
def callFimo(fastaFile, pwmFile, paraDict, outDir):
	"""Call the Fimo executable. Below are the parameters 
		USAGE: fimo [options] <motif file> <sequence file>
		The default command as run on the online site is:
		fimo --oc . --verbosity 1 --thresh 0.0001 pwm_File fasta_file

	   Options:
		 --alpha <double> (default 1.0)
		 --bgfile <background> (default from NR sequence database)
		 --max-seq-length <int> (default=2.5e8)
		 --max-stored-scores <int> (default=100000)
		 --max-strand
		 --motif <id> (default=all)
		 --motif-pseudo <float> (default=0.1)
		 --no-qvalue
		 --norc
		 --o <output dir> (default=fimo_out)
		 --oc <output dir> (default=fimo_out)
		 --parse-genomic-coord
		 --psp <PSP filename> (default none)
		 --prior-dist <PSP distribution filename> (default none)
		 --qv-thresh
		 --text
		 --thresh <float> (default = 1e-4)
		 --verbosity [1|2|3|4] (default 2)

	   Use '-' for <sequence file> to read the database from standard input.
	   Use '--bgfile motif-file' to read the background from the motif file.
		Args:
	
	Returns:
	
	"""
	#get the path of the module	
	Fimo_path = os.path.realpath(__file__)
	#remove the name of the module from the end of it. 
	split = Fimo_path.split('/')
	Fimo_path = '/'.join(split[:len(split)-1])
	Fimo_path = Fimo_path + '/fimo'
	command = (Fimo_path + ' ' + ' '.join(['--{} {}'.format(k,v) for k,v in paraDict.iteritems()]) + ' --oc ' + outDir + ' --verbosity 1' + 
		' ' + pwmFile + ' ' + fastaFile  )
	print '\n\nFIMO command:', command
	print '\nFIMO out folder:',outDir,'\n'
	try:
		os.system(command)
	except:
		print 'FIMO execution failed. Exiting'
		exit()



def parseFimo(fimoFile, strand):
	""" parse the fimo.txt file
	Args:
		the fimo.txt file
		strand = single or double
	Returns:
		fimoDict: a dict between motif ID and a list of sequences it occurs in 
	
	"""
	#dict to store for each motif list of seqs that it occurs in
	fimoDict = {}
	#read the fimo.txt file
	with open(fimoFile, 'rb') as handler:
		for line in handler:
			line = line.strip()
			#if re.search(r'#', line):
				#continue
			if re.search(r'stop', line):
				continue
			lineSplit = line.split()
			motifName = lineSplit[0]
			seqName = lineSplit[1]
			start = int(lineSplit[2])
			stop = int(lineSplit[3])
			pval = float(lineSplit[6])
			#check of this is on the negative strand hit
			if strand == 'single' and inStrand == '-':
				continue
			if motifName not in fimoDict:
				fimoDict[motifName] = []
			if seqName not in fimoDict[motifName]:
				fimoDict[motifName].append(seqName)
	
	#check the motifs
	print '\n\nfimo number of motifs:', len(fimoDict)
	#for motifName, seqList in fimoDict.items():
		#print motifName
		#print '\t', seqList
	
	return fimoDict
	

def main(args):
    pass

##
if( __name__ == '__main__' ):
    main(sys.argv)
#--eof--#

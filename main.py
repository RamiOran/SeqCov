# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
from __future__ import division
"""
main.py -- main program of motif sequence coverage  pipeline tool
example for running code: python main.py -jid test -confFile ./data/default_scan_only.conf
This will create a results folder with the name test. You will need to delete the folder if you want to run the code again with the same name
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
sys.path.append(inPath + '/algo')
import greedy
import motif_pwm_scan_cov
import motif_pwm_scan_only



#MAIN
def main(args):
	#example for running code: python main.py -jid test -confFile ./data/default_scan_only.conf
    print "main.py::main()"
    parser = argparse.ArgumentParser()
    parser.add_argument("-jid", "--jid",  help="enter job ID") #job id to make a folder to store all the data for a specific job
    parser.add_argument("-confFile", "--confFile",  help="enter the configuration file")#path to configuration file
    args = parser.parse_args()
    print 'jobId:', args.jid,'configfile:', args.confFile
    
	#make a results directory to store job results
    resultsDirName = args.jid
    os.makedirs(resultsDirName)
    
    #make a file list to store all the files to be moved to the results folder
    fileList = []
    
    #copy the config file
    cpConfFileName = args.jid + '_in_conf_file'
    cpConfFile = open(cpConfFileName, 'wb')
    with open(args.confFile, 'rb') as handler:
            for line in handler:
                cpConfFile.write(line)
    cpConfFile.close()
    fileList.append(cpConfFileName)		
	
    #make a config object
    confObj = conf.Conf()
    confDict = confObj.read(args.confFile)
    
    
    ############
    #Have a PWM file and want to scan it across a fasta file and then apply sequence coverage
    ############
    if confDict['job.type']['type'] == 'motifPwmScanCov':
        print 'motif_pwm scanning and coverage operation'
        motif_pwm_scan_cov.callMotifPwmScanCov(args, confDict, fileList)
        #move files to results folder
        for outFile in fileList:
            shutil.move(outFile, resultsDirName)
        exit()
     
     
    ############
    #Have a PWM file and want to scan it across a file only
    ############  
    if confDict['job.type']['type'] == 'pwmScan':
        print 'motif pwm scanning only'
        motif_pwm_scan_only.callMotifPwmScanOnly(args, confDict, fileList)
        #move files to results folder
        for outFile in fileList:
            shutil.move(outFile, resultsDirName)
        exit()
        
    ###############
    #EXIT
    ###############
    exit()
			
	
    

#calling main
if( __name__ == "__main__" ):
    main(sys.argv)
  

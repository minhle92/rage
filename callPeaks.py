import numpy as np
from Bio import SeqIO
import copy
import re

   
def parseGenome(genomeFile):
	"""Parses genome as a string and returns it along with its length"""
	print genomeFile
	handle=open(genomeFile)
	record = SeqIO.read(genomeFile, "fasta")
	return len(record.seq), str(record.seq)

def createHist(coords, genomeSize):
	"""creates a list storing the number of reads mapping to each base position
	in the genome"""
	print "creating hist"
	readCounts=[0]*genomeSize
	for i in xrange(len(coords)):
		start=coords[i][0]
		end=coords[i][1]
		for j in xrange(start, end+1):
			readCounts[j]+=1
	return readCounts


def readCoords(filename):
	"""Reads start and stop position of each read's mapping to the genome
	and parses it into a list for further processing"""
	with open(filename, "r") as f:
		contents = f.read()
	lines = contents.split("\n")
	coordList=[]
	for i in xrange(len(lines)-1):
		cols = lines[i].split()
		startIndex = "0"
		read = "0"
		for j in range(len(cols)): 
			if (j == 3):
				startIndex = cols[3]
			if (j == 9):
				read = cols[9]
		if read == "0": 
			continue
		else:
			if len(read) != 28:
				continue
		coordList.append([int(startIndex),int(startIndex) + len(read)])
	return coordList

def deleteZeros (input) :
	"""Removes peaks with zero values from the input"""
	result = []
	for i in xrange(len(input)): 
		if input[i] != 0:
			result.append(input[i])
	return result	

def getPeakThreshold(histogram) :
    
    nonzero_data = deleteZeros(histogram)

    median = np.median(nonzero_data)
    mean = np.mean(nonzero_data)

    distance = [0] * len(histogram)

    for i in xrange(len(histogram)):
        if histogram[i] != 0: 
            distance[i] = np.abs(histogram[i] - median)

    nonzero_distance = deleteZeros(distance)
        
    d_median = np.median(nonzero_distance)

    print "median: \n"
    print median
    # print "distance array: \n"
    # print distance
    print "distance median: \n"
    print d_median

    medianRatio = [0] * len(histogram)

    if d_median != 0:
        for i in xrange(len(histogram)):
            medianRatio[i] = distance[i]/d_median

    # print "medianRatio: "
    # print medianRatio

    outlierIndices = []
    #delete outliers from histogram array
    for i in xrange(len(histogram)):
        if medianRatio[i] > 2:
            histogram[i] = 0
            outlierIndices.append(i)

	#collect all non-zero data points after outliers have been zeroed to 
	#calculate standard deviation
    nonzero_noOutliers = deleteZeros(histogram)
    stddev = np.std(nonzero_noOutliers)

    print "standard deviation: \n"
    print stddev

	#once you have standard deviation, go back and fix outlier points to 
	#more reasonable values
    for i in xrange(len(outlierIndices)):
        histogram[outlierIndices[i]] = median + 2 * stddev
    
        # print "final histogram array: \n"
        # print histogram

    nonzero_fixedOutliers = deleteZeros(histogram)
    finalStddev = np.std(nonzero_fixedOutliers)
    finalMean = np.mean(nonzero_fixedOutliers)
    peakThresh = finalMean + 1 * finalStddev

    print "final Stddev: "
    print finalStddev
    print "final mean: "
    print finalMean
    print "peak threshold: "
    print peakThresh
    return peakThresh

def callPeaks(histogram, peakThresh, genomeSize):
	"""Generates a string of ones and zeroes to denote whether a given
	nucleotide position is considered a peak or not"""
	peakCalls=['0']*genomeSize
	for i in xrange(genomeSize):
		if histogram[i]>=peakThresh:
			peakCalls[i]='1'
	
	return ''.join(peakCalls)

def findPeriodicityR (inputStr, matchesList, minimumLength): 
    """Regular expression matching to find regions where there is periodicity
    of ribosome occupation peaks"""
    matches = list(re.finditer("((1..)+(...){0,10})+", inputStr))

    if len(matches) ==  0: 
        return matchesList
    else:
        inputList = list(inputStr)
        for m in matches: 
            start = m.start()
            end = m.end()
            for i in range(start, end, 3): 
                inputList[i] = "0"
            if (end - start) >= minimumLength:
                matchesList.append((start, end))

        return findPeriodicityR(''.join(inputList), matchesList, minimumLength)

def findPeriodicity (inputStr, minimumLength): 
    """Wrapper function"""
    return findPeriodicityR(inputStr, [], minimumLength)

	
def findCodonsInSegment(codons, start, stop, genomeString):
	"""function to find and return the first in frame codon in a specified codons list
	in the region defined by start and stop on the genome"""
	for i in xrange(start, stop, 3):
		if i+3<=len(genomeString):
			if genomeString[i:i+3] in codons:
				return i
			
	return None
	
def findOrfs(startStopList, genomeString, maxAllowance, startCodons, stopCodons):	
	"""Function to predict ORFs by looking for start and stop codons near blocks
	containing peak periodicity"""
	orfList=[]
	for (start, stop) in startStopList:
		orfStart=findCodonsInSegment(startCodons, start, stop, genomeString)
		if orfStart==None:		
			orfStart=findCodonsInSegment(startCodons, start-3*maxAllowance, start, genomeString)
		if orfStart!=None:
			orfEnd=findCodonsInSegment(stopCodons, orfStart, len(genomeString), genomeString)
			if orfEnd!=None:
				orfFrame=orfStart%3+1
				orfList.append([orfStart, orfEnd, orfFrame])
				
	return orfList

def createComplement(genomeString):
	"""Returns a complement of the genome string in proper orientation i.e. 5' to 3'"""
	complement={'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a':'t', 'c':'g', 'g':'c','t':'a'}
	#substitutions being made
	comp=[complement[base] for base in genomeString]
	compStr=''.join(comp)
	#return the reverse of it
	return compStr[::-1]

def myCMP(item1, item2):
	"""function to sort the orfList by start position of ORFs"""
	if item1[0]>item2[0]:
		return 1
	elif item1[0]<item2[0]:
		return -1
	else:
		return 0
	
	
def writeToFile(records, out_file, genomeSize):
	"""Writes a file similar to a .predict file for the predicted ORFs"""
	out_f=open(out_file, 'w')
	len_str=len(str(genomeSize))
	for i in xrange(len(records)):
		[start, stop, frame]=records[i]
		#output formatting and padding by required number of characters
		serialStr='{num:{fill}{width}}'.format(num=i+1, fill='0', width=len_str)
		startStr='{num:{fill}{width}}'.format(num=start+1, fill=' ', width=len_str)	
		stopStr='{num:{fill}{width}}'.format(num=stop, fill=' ', width=len_str)
		frameStr=frame
		out_f.write(serialStr+" "+startStr+" "+stopStr+" "+frameStr+"\n")
		
	out_f.close()

	
def run(in_file,genomeFile,out_file):
	#in_file="../../results/data/outputsfull/bowtie/ribotrimmedfiltered_rRNA_tRNA.sam" #input sam file
	#genomeFile="../../data/ecolik12/k12genome.fasta" #input genome fasta file
	#output .predict file
	#out_file="callPeaksTest-2.txt"
	#start codon list, can be specified by user as well
	startCodons=["ATG", "atg"]
	#stop codon list, can be specified by user as well
	stopCodons=["TAA", "taa", "TAG", "tag", "TGA", "tga"]
	#minimum peak periodicity gene length
	minBlockLength=300
	#maximum allowance before the block for the finding of a start codon
	maxCodonAllowance=50 #corresponds to 50*3 nucleotides
	genomeSize, genomeString=parseGenome(genomeFile)
	coords=readCoords(in_file)
	readCounts=createHist(coords, genomeSize)
	assert(genomeSize==len(readCounts))
	peakThresh=getPeakThreshold(readCounts)	
	peakString=callPeaks(readCounts, peakThresh, genomeSize)
	startStopPositive=findPeriodicity(peakString, minBlockLength)
	startStopNegative=findPeriodicity(peakString[::-1], minBlockLength)
	fwdOrfList=findOrfs(startStopPositive, genomeString, maxCodonAllowance, startCodons, stopCodons)
	revCompGenomeString=createComplement(genomeString)
	revOrfList=findOrfs(startStopNegative, revCompGenomeString, maxCodonAllowance, startCodons, stopCodons)
	for record in fwdOrfList:
		#proper formatting for writing to file
		record[2]="+"+str(record[2])
	for record in revOrfList:
		#array indices need to be subtracted by genomeSize on complement
		record[0]=genomeSize-1-record[0]
		record[1]=genomeSize-1-record[1]
		record[2]="-"+str(record[2])	
	#merge both lists
	fwdOrfList.extend(revOrfList)
	fwdOrfList.extend(revOrfList)
	orfList=fwdOrfList
	#sort both lists with own compare function, sorts by start position
	orfList=sorted(orfList, cmp=myCMP)
	#writes to out_file
	writeToFile(orfList, out_file, genomeSize)
	
#run()

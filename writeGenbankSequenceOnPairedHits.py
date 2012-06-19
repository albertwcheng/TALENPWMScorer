from sys import *
from BedSeqUtil import BedSeqClient
from GenbankFileWriter import GenbankFileWriter

def printUsageAndExit(programName):
	print >> stderr,"Usage: %s seqdir pairedhitFile flankingLength" %(programName)
	exit(1)


if __name__=='__main__':
	programName=argv[0]
	args=argv[1:]
	try:
		seqdir,filename,flankingLen=args
		flankingLen=int(flankingLen)
	except:
		printUsageAndExit(programName)
		
	
	bseq=BedSeqClient(seqdir,"bed")
	
	fil=open(filename)
	lino=0
	for lin in fil:
		lino+=1
		lin=lin.rstrip("\r\n")
		fields=lin.split("\t")
		pair1Name=fields[4]
		pair1Chrom=fields[1]
		pair1Start0=int(fields[2])
		pair1End1=int(fields[3])
		
		pair1Strand=fields[6]
		pair2Name=fields[11]
		pair2Start0=int(fields[9])
		pair2End1=int(fields[10])
		pair2Strand=fields[13]
		minStart0=min(pair1Start0,pair2Start0)
		maxEnd1=max(pair1End1,pair2End1)
		
		#print >> stderr,pair1Start0,pair1End1,pair2Start0,pair2End1,#flankingLen+pair1End1-pair1Start0-minStart0,pair1Strand
		flankingStart0=minStart0-flankingLen
		flankingEnd1=maxEnd1+flankingLen
		
		matchName=(pair1Name+"_"+pair2Name).replace(":","_")
		gbWriter=GenbankFileWriter(matchName)
		gbWriter.source=pair1Chrom+":"+str(minStart0+1)+"-"+str(maxEnd1)
		gbWriter.organism=pair1Chrom+":"+str(minStart0+1)+"-"+str(maxEnd1)
		gbWriter.sequence=bseq.getBedSeq(pair1Chrom+"\t"+str(flankingStart0)+"\t"+str(flankingEnd1)).split("\t")[3]
		#print >> stderr,gbWriter.sequence
		gbWriter.addFeature(flankingLen+pair1Start0-minStart0+1,flankingLen+pair1End1-minStart0,pair1Strand,"pair1")
		gbWriter.addFeature(flankingLen+pair2Start0-minStart0+1,flankingLen+pair2End1-minStart0,pair2Strand,"pair2")
		gbWriter.writeGenbankFile(str(lino)+"_"+matchName+".gb")
	
	fil.close()
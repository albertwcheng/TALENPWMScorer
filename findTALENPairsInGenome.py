#!/usr/bin/env python

from sys import *
from os import listdir,makedirs,system
from os.path import basename,exists,dirname

def printUsageAndExit(programName):
	print >> stderr,"Usage: %s seqDir mode1 motifDef1 modeParam1 topN1 minSpacer maxSpacer mode2 motifDef2 modeParam2 topN2 outDir" %(programName)
	exit(1)
	
def oppostrand(s):
	if s=="+":
		return "-"
	elif s=="-":
		return "+"
	
if __name__=='__main__':
	programName=argv[0]
	args=argv[1:]
	programDir=dirname(programName)
	TALENPWMScorerName="TALENPWMScorer"
	if programDir=="":
		programDir="."
	TALENPWMScorerPath=programDir+"/"+TALENPWMScorerName
	try:
		seqDir,mode1,motifDef1,modeParam1,topN1,minSpacer,maxSpacer,mode2,motifDef2,modeParam2,topN2,outDir=args
		if mode1 not in ["consensus","PWM"]:
			print >> stderr,"mode1",mode1,"not defined"
			raise Exception("")
		if mode2 not in ["consensus","PWM"]:
			print >> stderr,"mode2",mode2,"not defined"
			raise Exception("")
		modeParam1=int(modeParam1)
		modeParam2=int(modeParam2)
		topN1=int(topN1)
		topN2=int(topN2)
		minSpacer=int(minSpacer)
		maxSpacer=int(maxSpacer)
	except:
		printUsageAndExit(programName)
		
	
	if mode1=="consensus":
		w1=len(motifDef1)
	
	if mode2=="consensus":
		w2=len(motifDef2)
	
	filesInSeqDir=listdir(seqDir)
	
	refs=[]
	
	if not exists(outDir):
		makedirs(outDir)
	
	fchrom1filename=outDir+"/"+"chroms1.txt"
	fchroms=open(fchrom1filename,"w")
	
	
	
	for filename in filesInSeqDir:
		ref=basename(filename).replace(".seq","")
		refs.append(refs)
		print >> fchroms,ref
	
	fchroms.close()
	
	TALEN1foundOut=outDir+"/TALEN1.found"
	
	TALEN1Stderr=outDir+"/TALEN1.stderr"
	#now run the first TALEN
	command="%s %s %s %s %s %d %d > %s 2> %s" %(TALENPWMScorerPath,seqDir,fchrom1filename,mode1,motifDef1,modeParam1,topN1,TALEN1foundOut,TALEN1Stderr)
	
	print >> stderr,command
	
	system(command)
	

	TALEN1foundOutWithChrom2=TALEN1foundOut+".wChroms2"
	#now form the chroms2.txt
	fTALEN1found=open(TALEN1foundOut)
	fTALEN1foundWithChrom2=open(TALEN1foundOutWithChrom2,"w")
	fchrom2filename=outDir+"/chroms2.txt"
	fchroms=open(fchrom2filename,"w")
	for lin in fTALEN1found:
		lin=lin.rstrip("\r\n")
		if len(lin)==0 or lin[0]=="#":
			continue
		fields=lin.split("\t")
		ref=fields[1]
		gstart10=int(fields[2])
		gend11=int(fields[3])
		strand1=fields[6]
		strand2=oppostrand(strand1)
		if strand1=="+":
			gstart20=gend11+minSpacer
			gend21=gend11+maxSpacer+w2
		elif strand1=="-":
			gstart20=gstart10-maxSpacer-w2
			gend21=gstart10-minSpacer
			
		#chroms are start1 end1
		gstart21=gstart20+1
		confine="%s:%d-%d:%s" %(ref,gstart21,gend21,strand2)
		print >> fchroms,confine
		print >> fTALEN1foundWithChrom2,lin+"\t"+confine
		
	fchroms.close()
	fTALEN1found.close()
	fTALEN1foundWithChrom2.close()
	
	TALEN2foundOut=outDir+"/TALEN2.found"
	TALEN2Stderr=outDir+"/TALEN2.stderr"
	command="%s %s %s %s %s %d %d > %s 2> %s" %(TALENPWMScorerPath,seqDir,fchrom2filename,mode2,motifDef2,modeParam2,topN2,TALEN2foundOut,TALEN2Stderr)
	
	print >> stderr,command
	system(command)
	
	PairedHitOut=outDir+"/paired_hits.txt"
	PairedHitStderr=outDir+"/paired_hits.stderr"
	command="joinu.py -1 8 -2 1 %s %s > %s 2> %s" %(TALEN1foundOutWithChrom2,TALEN2foundOut,PairedHitOut,PairedHitStderr)
	print >> stderr,command
	system(command)
	
	
	PairedHitOutSimp=outDir+"/paired_hits_simp.txt"
	command="cuta.py -f5,12 %s > %s" %(PairedHitOut,PairedHitOutSimp)
	print >> stderr,command
	system(command)
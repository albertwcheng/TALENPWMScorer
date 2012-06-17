#!/usr/bin/env python

from sys import *
from optparse import OptionParser
from math import log
####diresidue counts were adapted from Science paper http://www.sciencemag.org/content/suppl/2009/10/29/1178817.DC1/Moscou.SOM.rev1.pdf

diresidue_counts={}
unknown_diresidue=[1,1,1,1]

def initModel(options):
	global diresidue_counts
	
	if options.model == 2:
		# known
		diresidue_counts['HD'] = [ 6, 45,  0,  1]
		diresidue_counts['NI'] = [35,  4,  0,  0]
		diresidue_counts['NG'] = [ 3,  4,  0,  32]
		diresidue_counts['NN'] = [ 7,  3,  8,  1]
		diresidue_counts['NS'] = [11,  3,  2,  0]
		diresidue_counts['N*'] = [ 1,  7,  0,  1]
		diresidue_counts['HG'] = [ 0,  2,  0,  3]
		diresidue_counts['IG'] = [ 0,  0,  0,  1]
		diresidue_counts['H*'] = [ 0,  0,  0, 1]
		# unknown
		#diresidue_counts['ND'] = [1, 1, 1, 1]
		#diresidue_counts['HA'] = [1, 1, 1, 1]
		#diresidue_counts['HI'] = [1, 1, 1, 1]
		#diresidue_counts['HN'] = [1, 1, 1, 1]
		#diresidue_counts['NA'] = [1, 1, 1, 1]
		#diresidue_counts['NK'] = [1, 1, 1, 1]
		#diresidue_counts['S*'] = [1, 1, 1, 1]
		#diresidue_counts['NH'] = [1, 1, 1, 1]
		#diresidue_counts['NC'] = [1, 1, 1, 1]
		#diresidue_counts['YG'] = [1, 1, 1, 1]
		#diresidue_counts['HH'] = [1, 1, 1, 1]
		#diresidue_counts['SN'] = [1, 1, 1, 1]
		#diresidue_counts['SS'] = [1, 1, 1, 1]
	elif options.model == 3:
		# known
		diresidue_counts['HD'] = [  7,  99,   0,   1]
		diresidue_counts['NI'] = [ 58,   6,   0,   0]
		diresidue_counts['NG'] = [  6,   6,   1,  57]
		diresidue_counts['NN'] = [ 21,   8,  26,   2]
		diresidue_counts['NS'] = [ 20,   6,   4,   0]
		diresidue_counts['N*'] = [  1,  11,   1,   7]
		diresidue_counts['HG'] = [  1,   2,   0,  15]
		diresidue_counts['HA'] = [  1,   4,   1,   0]
		diresidue_counts['ND'] = [  0,   4,   0,   0]
		diresidue_counts['NK'] = [  0,   0,   2,   0]
		diresidue_counts['HI'] = [  0,   1,   0,   0]
		diresidue_counts['HN'] = [  0,   0,   1,   0]
		diresidue_counts['NA'] = [  0,   0,   1,   0]
		diresidue_counts['IG'] = [  0,   0,   0,   1]
		diresidue_counts['H*'] = [ 0,  0,  0, 1]
		# unknown
		#diresidue_counts['S*'] = [1, 1, 1, 1]
		#diresidue_counts['NH'] = [1, 1, 1, 1]
		#diresidue_counts['YG'] = [1, 1, 1, 1]
		#diresidue_counts['SN'] = [1, 1, 1, 1]
		#diresidue_counts['SS'] = [1, 1, 1, 1]
		#diresidue_counts['NC'] = [1, 1, 1, 1]
		#diresidue_counts['HH'] = [1, 1, 1, 1]
	
	#add pseudocounts
	if options.pseudoCounts>0:
		for diresidue,counts in diresidue_counts.items():
			for i in range(0,4):
				counts[i]+=options.pseudoCounts
	
		
	if options.outMode=="prob" or options.outMode=="logprob" or options.outMode=="logkp":
		for diresidue,counts in diresidue_counts.items():
			countSum=sum(counts)
			for i in range(0,4):
				counts[i]=float(counts[i])/countSum*options.obsWeight+0.25*(1-options.obsWeight)
				if options.outMode=="logprob":
					counts[i]=log(counts[i],options.logb)
				elif options.outMode=="logkp":
					counts[i]=int(log(counts[i],options.logb)*1000)
		
		
	
		
	
#### END 

def printPWMOnDiresidues(stream,diresidues):

	consensus=""
	
	global diresidue_counts,unknown_diresidue
	#PWM=[]
	for diresidue in diresidues:
		if diresidue in diresidue_counts:
			#PWM.append(diresidue_counts[diresidue])
			prob=diresidue_counts[diresidue]
		else:
			#PWM.append(unknown_diresidue)
			prob=unknown_diresidue
		
		print >> stream,str(prob[0])+"\t"+str(prob[1])+"\t"+str(prob[2])+"\t"+str(prob[3])
		
		Argmax=-1
		if prob[0]==prob[1] and prob[1]==prob[2] and prob[2]==prob[3]:
			consensus+="N"
		else:
			maxprob=max(prob)
			if prob[0]==maxprob:
				consensus+="A"
			elif prob[1]==maxprob:
				consensus+="C"
			elif prob[2]==maxprob:
				consensus+="G"
			elif prob[3]==maxprob:
				consensus+="T"
		
	print >> stream,"#!consensus "+consensus
		
if __name__=='__main__':
	parser=OptionParser(usage="Usage: %prog [options] filename")
	parser.add_option("--model",dest="model",default=3,type="int",help="set model [2,*3]")
	parser.add_option("--out-mode",dest="outMode",default="logkp",type="string",help="set output mode, count: counts, prob: probabilities, logprob: log probabilities, *logkp: log probability score times 1000")
	parser.add_option("--out-dir",dest="outDir",default=".",type="string",help="set out dir default: current dir[.]")
	parser.add_option("--psuedo-count",dest="pseudoCounts",default=0,type="int",help="set pseudocount to be added to each base for each column [default:0]")
	parser.add_option("--obs-weight",dest="obsWeight",default=0.99,type="float",help="set the weight given to observed probabilities. finalPWM=PWM*w+0.25*(1-w) for each cell [default=0.99]")
	parser.add_option("--log-base",dest="logb",default=2,type="int",help="set the log base for logprob [default: 2]")

	(options,args)=parser.parse_args()
	
	try:
		filename,=args
	except:
		parser.print_help()
		exit(1)

	if options.outMode not in ["count","prob","logprob","logkp"]:
		print >> stderr,"--out-mode must be either count , prob, logprob"
		parser.print_help()
		exit(1)
			
	initModel(options)
	

	
	fil=open(filename)
	for lin in fil:
		if len(lin)<1 or lin[0]=="#":
			continue
		lin=lin.rstrip("\r\n")
		fields=lin.split()
		if len(fields)<2:
			continue
			
		TALENName=fields[0]
		diresidues=fields[1:]
		
		fout=open(options.outDir+"/"+TALENName+".rowmat","w")
		print >> fout,"#"+lin+" model="+str(options.model)+" outMode="+options.outMode+" pseudocount="+str(options.pseudoCounts)+(" obsweight="+str(options.obsWeight)+" logbase="+str(options.logb) if options.outMode!="count" else "")
		
		printPWMOnDiresidues(fout,diresidues)
		fout.close()

	fil.close()
	
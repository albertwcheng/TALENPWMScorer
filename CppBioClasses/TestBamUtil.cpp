#include <iostream>
#include <StringUtil.h>
#include "BamUtil.h"
using namespace std;

int main(int argc,char **argv)
{
	if(argc<5){
		cerr<<argv[0]<<" "<<"bamfile chrom start end"<<endl;
		return 1;	
	}
	
	
	//int num=BamReader::countTotalNumOfReadsInBam(argv[1]);
	//cerr<<"number of count in "<<argv[1]<<" is "<<num<<endl;
	
	BamReader bamReader(argv[1]);
	string chrom=argv[2];
	int start1=StringUtil::atoi(argv[3]);
	int end1=StringUtil::atoi(argv[4]);
	string genomBCoord=BamReader::getGenomeBrowserCoordinate(chrom,start1,end1);
	int numReads=bamReader.fetchCountOverlappingRegion(chrom, start1-1,  end1,true);
	bamReader.close();
	
	cerr<<"number of reads in "<<genomBCoord<<" is "<<numReads<<endl;
	
	return 0;	
}
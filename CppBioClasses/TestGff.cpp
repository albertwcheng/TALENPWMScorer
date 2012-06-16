#include <iostream>
#include "Gff.h"
#include <set>

using namespace std;
using namespace Gff;


int main(int argc,char **argv)
{
	/*set<SmartPtr<Gene> > geneSet;
	geneSet.insert(new Gene("chr1",0,100,"geneB"));
	geneSet.insert(new Gene("chr1",0,100,"geneA"));
	geneSet.insert(new Gene("chr1",0,100,"geneC"));
	geneSet.insert(new Gene("chr1",50,100,"geneZ"));
	geneSet.insert(new Gene("chr1",50,60,"geneY"));
	
	for(set<SmartPtr<Gene> >::iterator i=geneSet.begin();i!=geneSet.end();i++){
		cerr<<(*i)->chrom()<<"\t"<<(*i)->start0()<<"\t"<<(*i)->end1()<<"\t"<<(*i)->name()<<endl;
		delete *i;
	}*/
	
	if(argc<2){
		cerr<<argv[0]<<"bedfile"<<endl;
		return 1;
	}
	
	
	Annotation annot;
	annot.readBedFile(argv[1]);
	annot.print(cerr);
	
	annot.name_genes["ENAH"]->constructBlockInfo();
	annot.name_genes["ENAH"]->printBlockInfo(cerr);
	
	
	//annot.writeBed(cout);
	return 0;
}
#include <string>
#include <set>

#include <iostream>

using namespace std;

#include <StringUtil.h>
#include <SortPair.h>
#include <Array.h>
#include <limits.h>
#include <SmartPtr.h>
#include <fstream>

#define STRAND_FORWARD '+'
#define STRAND_REVERSE '-'
#define STRAND_UNKNOWN '.'


#define INV_START0 INT_MAX
#define INV_END1 INT_MIN


#define MIN(x,y) ((x)<(y)?(x):(y))
#define MAX(x,y) ((x)>(y)?(x):(y))


namespace Gff {


	typedef int Coord;
	
	class Bound01:public SortPair<Coord,Coord>
	{
		// k1: start0
		// k2: end1
		public:
		syn_rwaccessors(start0,k1,Coord);
		syn_rwaccessors(end1,k2,Coord);
		inline Bound01(int _start0=INV_START0,int _end1=INV_END1){
			set_start0(_start0);
			set_end1(_end1);
		}
	};
	
	class ChrBound01:public SortPair<string,Bound01>
	{
		//k1: chrom
		//k2: bound01
		// k2.k1: start0
		// k2.k2: end1
		
		public:
		syn_rwaccessors(chrom,k1,string);
		syn_rwaccessors(bound,k2,Bound01);
		
		syn_rwaccessors(start0,k2.k1,Coord);
		syn_rwaccessors(end1,k2.k2,Coord);
		
		inline ChrBound01(const string& _chrom="",Coord _start0=INV_START0,Coord _end1=INV_END1){
			set_chrom(_chrom);
			set_start0(_start0);
			set_end1(_end1);
		}
		
		inline ChrBound01(const string& _chrom,const Bound01& _bound){
			set_chrom(_chrom);
			set_bound(_bound);
		}
		
	};
	
	
	class Transcript;
	
	class TranscriptExonRank
	{
	public:
		SmartPtr<Transcript> parent;
		int rank;
		inline TranscriptExonRank(SmartPtr<Transcript> _parent,int _rank):parent(_parent),rank(_rank){
			
		}
	};
	
	class Exon:public Bound01
	{
		
	public:
		//TODO
		Array<string>* names;
		Array<string>* IDs;
		
		syn_rwaccessors(start0,k1,Coord);
		syn_rwaccessors(end1,k2,Coord);
		
		
		Array<TranscriptExonRank>* transcriptLinks;
		inline int length(){
			return end1()-start0();
		}
		inline Exon(const Bound01 &_bound):Bound01(_bound),names(NULL),IDs(NULL),transcriptLinks(NULL)
		{
			
		}
		
		
		
	
		inline string getFirstName(){
			if(!names || names->size()==0){
				return "";
			}
			return (*names)[0];
		}
		inline  string getFirstID(){
			if(!IDs || IDs->size()==0){
				return "";
			}
			return (*IDs)[0];
		}
		
		inline void addNameAndID(const string& _name,const string& _id){
			if(!names){
				names=new Array<string>;
			
			}
			if(!IDs){
				IDs=new Array<string>;
			}
			
			names->push_back(_name);
			IDs->push_back(_id);
		}
		
		inline ~Exon(){
			if(names){
				delete names;
			}
			if(IDs){
				delete IDs;
			}
			if(transcriptLinks){
				delete transcriptLinks;
			}
		}
		
		inline void print(ostream& os,const string&prefix="");
				
	};
	
	
	class Transcript:public Bound01
	{
		//transcripts don't own exons. exons are owned by genes.
		public:
		
		string name;
		string ID;
		int thickStart;
		int thickEnd;
		string itemRgb;
		string score;
		
		syn_rwaccessors(start0,k1,Coord);
		syn_rwaccessors(end1,k2,Coord);
		
		set<SmartPtr<Exon> > exons;
		
		inline void addExon(SmartPtr<Exon> _exon,int _rank){
			exons.insert(_exon);
			//now, register the rank struct in exon
			if(!_exon->transcriptLinks){
				_exon->transcriptLinks=new Array<TranscriptExonRank>;
			}
			
			_exon->transcriptLinks->push_back(TranscriptExonRank(this,_rank));
			//now expand my bound
			set_start0(MIN(_exon->start0(),start0()));
			set_end1(MAX(_exon->end1(),end1()));
			
		}
		
		inline void print(ostream& os,const string& prefix=""){
			os<<prefix<<(start0()+1)<<"\t"<<end1()<<"\t"<<name<<"\t"<<ID<<endl;
			for(set<SmartPtr<Exon> >::iterator i=exons.begin();i!=exons.end();i++){
				(*i)->print(os,prefix+"\t");
			}
		}
		
	};

	
	class BlockInfo{
	public:
		Array<Coord> blockBounds;
		Array<Coord> blockCounts;
	};
	
	class Gene:public SortPair<ChrBound01,string>
	{
		//k1: ChrBound01
		// k1.k1: chrom
		// k1.k2: Bound01
		// k1.k2.k1: start0
		// k1.k2.k2: end1
		//k2: name
		public:
		
		syn_rwaccessors(chrbound,k1,ChrBound01);
		syn_rwaccessors(name,k2,string);
		syn_rwaccessors(chrom,k1.k1,string);
		syn_rwaccessors(bound,k1.k2,Bound01);
		syn_rwaccessors(start0,k1.k2.k1,Coord);
		syn_rwaccessors(end1,k1.k2.k2,Coord);
		
						
		Array<SmartPtr<Transcript> > transcripts;
		set<SmartPtr<Exon> > exons;
		char strand;
		
		BlockInfo blockInfo; 
		
		inline void constructBlockInfo(bool forceReconstruct=false)
		{
			if(blockInfo.blockBounds.size()>0 && !forceReconstruct){
				//already constructed
				return;
			}
			
			map<Coord,int> blockBounds;
			for(set<SmartPtr<Exon> >::iterator i=exons.begin();i!=exons.end();i++){
				blockBounds.insert(map<Coord,int>::value_type((*i)->start0(),0));
				blockBounds.insert(map<Coord,int>::value_type((*i)->end1(),0));
			}
			
			
			for(set<SmartPtr<Exon> >::iterator i=exons.begin();i!=exons.end();i++){
				map<Coord,int>::iterator startB=blockBounds.find((*i)->start0());
				map<Coord,int>::iterator endB=blockBounds.find((*i)->end1());
				for(map<Coord,int>::iterator j=startB;j!=endB;j++){
					j->second+=(*i)->transcriptLinks->size();
				}
			}
			
			//now convert to more convenient struct
			blockInfo.blockBounds.clear();
			blockInfo.blockCounts.clear();
			for(map<Coord,int>::iterator i=blockBounds.begin();i!=blockBounds.end();i++){
				blockInfo.blockBounds.push_back(i->first);
				blockInfo.blockCounts.push_back(i->second);
			}
		}
		
		inline void printBlockInfo(ostream& os){
			for(int i=0;i<blockInfo.blockBounds.size()-1;i++){
				int thisBlockStart0=blockInfo.blockBounds[i];
				int thisBlockEnd1=blockInfo.blockBounds[i+1];
				int thisBlockCount=blockInfo.blockCounts[i];
				os<<chrom()<<"\t"<<(thisBlockStart0)<<"\t"<<thisBlockEnd1<<"\t"<<thisBlockCount<<endl;
			}
		}
		
		
		
		inline Gene(const ChrBound01& _chrbound01, const string& _name="", char _strand=STRAND_UNKNOWN){
			set_chrbound(_chrbound01);
			set_name(_name);
			strand=_strand;
		}
		
		inline Gene(const string& _chrom="", Coord _start0=INV_START0, Coord _end1=INV_END1, const string& _name="", char _strand=STRAND_UNKNOWN){
			set_chrom(_chrom);
			set_start0(_start0);
			set_end1(_end1);
			set_name(_name);
			strand=_strand;
		}
		
		
		inline void addTranscript(SmartPtr<Transcript> _transcript){
			transcripts.push_back(_transcript);
			//expand gene bound?
			set_start0(MIN(_transcript->start0(),start0()));
			set_end1(MAX(_transcript->end1(),end1()));
			
		}
		
		inline SmartPtr<Exon> makeExonAvailable(int _start0,int _end1){
			SmartPtr<Exon> newExon=new Exon(Bound01(_start0,_end1));
			pair<set<SmartPtr<Exon> >::iterator,bool> result=exons.insert(newExon);
			if(!result.second){
				//exon not inserted, existed: delete the temp
				delete newExon;
				
				return *result.first;
				
			}else{
				return newExon;
			}
			
			
		}
		
		
		inline ~Gene(){
			//cerr<<"destroying gene"<<name()<<endl;
			for(Array<SmartPtr<Transcript> >::iterator i=transcripts.begin();i!=transcripts.end();i++){
				delete *i;
			}
			for(set<SmartPtr<Exon> >::iterator i=exons.begin();i!=exons.end();i++){
				delete *i;
			}
		}
		
		inline void print(ostream& os,const string& prefix=""){
			os<<prefix<<chrom()<<"\t"<<(start0()+1)<<"\t"<<end1()<<"\t"<<strand<<"\t"<<name()<<endl;
			for(Array<SmartPtr<Transcript> >::iterator i=transcripts.begin();i!=transcripts.end();i++){
				(*i)->print(os,prefix+"\t");
			}
		}
		
	};
	
	class Annotation{
	public:
		Array<SmartPtr<Gene> > genes;
		map<string,SmartPtr<Gene> > name_genes;
		map<string/*chrom*/, set<SmartPtr<Gene> >* > sorted_genes;
		typedef map<string,SmartPtr<Gene> >::iterator NameGeneMapI;
		
		inline NameGeneMapI name_genes_begin(){
			return name_genes.begin();
		}
		
		inline NameGeneMapI name_genes_end(){
			return name_genes.end();
		}
		
		inline string getGeneNameSepByDot(string transcriptName){
			vector<string> nameComp;
			StringUtil::split(transcriptName,".",nameComp);
			if(nameComp.size()==1){
				return transcriptName;
			}
			
			nameComp.resize(nameComp.size()-1); //is there a more efficient way?
			return StringUtil::join<vector<string>,vector<string>::const_iterator>(nameComp,".");
		}
		
		inline ~Annotation()
		{
			for(Array<SmartPtr<Gene> >::iterator i=genes.begin();i!=genes.end();i++){
				delete *i;
			}
			
			for(map<string/*chrom*/, set<SmartPtr<Gene> >* >::iterator iSortedGeneChrom=sorted_genes.begin();iSortedGeneChrom!=sorted_genes.end();iSortedGeneChrom++){
			
				delete iSortedGeneChrom->second;
			}
		}	
		
		inline void print(ostream& os,const string& prefix=""){
			for(map<string,SmartPtr<Gene> >::iterator i=name_genes.begin();i!=name_genes.end();i++){
				i->second->print(os,prefix);
			}
		}
		
		
		inline void readBedFile(string filename,string fieldSeparator="\t"){
			ifstream fil(filename.c_str());
			
			vector<string> fields;
			vector<int> blockSizes;
			vector<int> blockStarts;
			
			string prevChrom="";
			set<SmartPtr<Gene> >* chrGeneSet=NULL;
			
			while(fil.good()){
				string line;
				getline(fil,line);
				//1: chrom
				//2: chromStart
				//3: chromEnd
				//4: name
				//5: score
				//6: strand
				//7: thickStart
				//8: thickEnd
				//9: itemRgb
				//10: blockCount
				//11: blockSizes
				//12: blockStarts
				if(line.length()<0 || line[0]=='#'){
					continue;
				}
				
				
				
				StringUtil::split(line,fieldSeparator,fields);
				
				if(fields.size()<12){
					continue;
				}
				
				string& chrom=fields[0];
				int start0=StringUtil::atoi(fields[1]);
				int end1=StringUtil::atoi(fields[2]);
				string& transcriptName=fields[3];
				string& ID=fields[3];
				
				string geneName=getGeneNameSepByDot(transcriptName);
				
				string& score=fields[4];
				char strand=fields[5][0];
				int thickStart=StringUtil::atoi(fields[6]);
				int thickEnd=StringUtil::atoi(fields[7]);
				string& itemRgb=fields[8];
				int blockCount=StringUtil::atoi(fields[9]);
				StringUtil::splitInt2(fields[10],",",blockSizes);
				StringUtil::splitInt2(fields[11],",",blockStarts);
				
				
				//does gene exist already?
				SmartPtr<Gene> thisGene=NULL;
				
				
				
				map<string,SmartPtr<Gene> >::iterator iGene=name_genes.find(geneName);
				if(iGene==name_genes.end()){
					//no!
					thisGene=new Gene(chrom,start0,end1,geneName,strand);
					name_genes.insert(map<string,SmartPtr<Gene> >::value_type(geneName,thisGene));
					genes.push_back(thisGene);
					
					
					
					if(prevChrom!=chrom || chrGeneSet==NULL ){
						//not the same chrom as previous one request another 
						map<string/*chrom*/, set<SmartPtr<Gene> >* >::iterator iSortedGeneChrom=sorted_genes.find(chrom);
						if(iSortedGeneChrom==sorted_genes.end()){
							//need new one
							chrGeneSet=new set<SmartPtr<Gene> >;
							sorted_genes.insert(map<string/*chrom*/, set<SmartPtr<Gene> >* >::value_type(chrom,chrGeneSet));
							
						}else{
							chrGeneSet=iSortedGeneChrom->second;
						}
						
						prevChrom=chrom;
						
						
					}
					
					
					chrGeneSet->insert(thisGene);
					
					
					
				}else{
					thisGene=iGene->second;
				}
				
				SmartPtr<Transcript> newTranscript=new Transcript;
				newTranscript->name=transcriptName;
				newTranscript->ID=transcriptName;
				newTranscript->itemRgb=itemRgb;
				newTranscript->score=score;
				newTranscript->thickStart=thickStart;
				newTranscript->thickEnd=thickEnd;
				
				thisGene->addTranscript(newTranscript);
				
				for(int i=0;i<blockCount;i++){
					int exonStart0=start0+blockStarts[i];
					int exonEnd1=exonStart0+blockSizes[i];
					SmartPtr<Exon> thisExon=thisGene->makeExonAvailable(exonStart0,exonEnd1);
					//add exon to transcript
					newTranscript->addExon(thisExon,i);
				}
			}
			fil.close();
		}
		
		inline void writeBed(ostream& os){
			for(map<string/*chrom*/, set<SmartPtr<Gene> >* >::iterator i=sorted_genes.begin();i!=sorted_genes.end();i++){
				set<SmartPtr<Gene> > * chrGeneSet=i->second;
				for(set<SmartPtr<Gene> >::iterator j=chrGeneSet->begin();j!=chrGeneSet->end();j++){
					SmartPtr<Gene> thisGene=*j;
					for(Array<SmartPtr<Transcript> >::iterator k=thisGene->transcripts.begin();k!=thisGene->transcripts.end();k++){
						SmartPtr<Transcript> transcript=*k;
						os<<thisGene->chrom()<<"\t"<<transcript->start0()<<"\t"<<transcript->end1()<<"\t"<<transcript->name<<"\t"<<transcript->score<<"\t"<<thisGene->strand<<"\t"<<transcript->thickStart<<"\t"<<transcript->thickEnd<<"\t"<<transcript->itemRgb<<"\t";
						os<<transcript->exons.size()<<"\t";
						
						vector<string> blockSizes;
						vector<string> blockStarts;
						
						for(set<SmartPtr<Exon> >::iterator l=transcript->exons.begin();l!=transcript->exons.end();l++){
							SmartPtr<Exon> thisExon=*l;
							blockSizes.push_back(StringUtil::str(thisExon->length()));
							blockStarts.push_back(StringUtil::str(thisExon->start0()-transcript->start0()));
						}
						
						os<<StringUtil::join<vector<string>,vector<string>::const_iterator>(blockSizes,",")<<","<<"\t";
						os<<StringUtil::join<vector<string>,vector<string>::const_iterator>(blockStarts,",")<<","<<endl;
						
						
					}
				}
			}
		}
	};
	
	
	inline void Exon::print(ostream& os,const string&prefix){
		os<<prefix<<(start0()+1)<<"\t"<<end1();
		if(transcriptLinks){
			os<<"\t"<<"[";
			os<<transcriptLinks->at(0).parent->name<<"@"<<transcriptLinks->at(0).rank;
			for(int i=1;i<transcriptLinks->size();i++){
				os<<","<<transcriptLinks->at(i).parent->name<<"@"<<transcriptLinks->at(i).rank;
			}
		}
		os<<"]"<<endl;
	}

};	

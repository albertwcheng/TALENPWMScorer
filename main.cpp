#include "CppUtilClasses/Matrix.h"
#include "CppUtilClasses/StringUtil.h"
#include <iostream>
#include <string>
#include <fstream>
#include <limits>
#include <cctype>
#include <map>
#include "bedSeq/RandomAccessFile.h"
#include "CppBioClasses/Nucleic.h"


#define FORWARD_STRAND '+'
#define REVERSE_STRAND '-'
#define BOTH_STRANDS '.'

#define PWMFINDER_BUFFER_SIZE 1024

using namespace std;

class PosMatrix
{
	public:
	Matrix<double> matData;
	void readFile(const string&filename){
			
	}
};


class KmerScorer{
	public:
	virtual int score(const char* _kmer) const{return 0;}
	virtual int getLength() const{return 0;}
	virtual string toString() const{return "";}
};

class SimpleConsensusKmerScorer:public KmerScorer{
	private:
	int mismatch;
	string consensus;
	public:
	SimpleConsensusKmerScorer(const string& _consensus, int _mismatch):consensus(_consensus),mismatch(_mismatch){}
	string toString() const{
		return consensus+"/"+StringUtil::str(mismatch);	
	}
	int score(const char*_kmer) const{
		const char*consensus_c_str=consensus.c_str();
		int consensus_length=consensus.length();
		int mismatches=0;
		for(int i=0;i<consensus_length;i++){
			if (toupper(consensus[i])!=toupper(_kmer[i]))
			{
				mismatches++;
				if (mismatches>mismatch){
					return  numeric_limits<int>::min();
				}
			}
		}
		
		return int(float(consensus_length-mismatches)/consensus_length*10000);
	}
	
	int getLength() const
	{
		return consensus.length();
	}
	
	SimpleConsensusKmerScorer createReverseComplementScorer() const;
	
};


SimpleConsensusKmerScorer SimpleConsensusKmerScorer::createReverseComplementScorer() const
{
	return SimpleConsensusKmerScorer(reverse_complement(this->consensus),this->mismatch);	
}

class PWMHitsFinder
{
	
	public:
	
	class Hit
	{
		public:
		int score;
		string ref;
		int start0;
		int end1;
		char strand;
		Hit(int _score,const string& _ref,int _start0,int _end1,char _strand):score(_score),ref(_ref),start0(_start0),end1(_end1),strand(_strand){}
		void print(ostream& os) const{
			string name=ref+":"+StringUtil::str(start0+1)+"-"+StringUtil::str(end1)+":"+strand;
			os<<ref<<"\t"<<start0<<"\t"<<end1<<"\t"<<name<<"\t"<<score<<"\t"<<strand<<endl;
		}
	};
	
	int topN;
	KmerScorer* scorer;
	KmerScorer* reverseScorer;
	int bufferSize;
	multimap<int,Hit,greater<int> > hits;
	
	PWMHitsFinder(KmerScorer* _scorer,KmerScorer* _reverseScorer,int _topN,int _bufferSize):scorer(_scorer),reverseScorer(_reverseScorer),topN(_topN),bufferSize(_bufferSize){}
	
	void printHits(ostream& os){
		for(multimap<int,Hit,greater<int> >::iterator i=hits.begin();i!=hits.end();i++){
			(*i).second.print(os);
		}
	}
	
	
	void proposeHit(const Hit& hit)
	{
		if (hit.score==numeric_limits<int>::min())
			return;
			
		hits.insert(map<int,Hit,greater<int> >::value_type(hit.score,hit));
		if(topN>0){
			if(hits.size()>topN){
				//cerr<<"try erase"<<endl;
				hits.erase((++hits.rbegin()).base()); //remove the (topN+1)-ranked hit 			
				//cerr<<"erase ok"<<endl;	
			}
		}
		
	
	}
	
	void scanSeq(RandomAccessFile & raf,const string& ref,int gstart0,int gend1,char strand)
	
	
	{
		
		
		if(!scorer){
			if(strand==FORWARD_STRAND || strand==BOTH_STRANDS){
				cerr<<"forward scorer not specified"<<endl;
				return;
			}
		}	
		
		if(!reverseScorer){
			if(strand==REVERSE_STRAND || strand==BOTH_STRANDS){
				cerr<<"reverse scorer not specified"<<endl;
				return;	
			}	
		}
		
		
		
		string buffer;
		int curgstart0=gstart0;
		
		
		buffer=raf.get(curgstart0,min(curgstart0+bufferSize,gend1));

		const char*buffer_c_str=buffer.c_str();
		int k=scorer->getLength();
		int fedSize=buffer.length();
		const char*curKmer;


		
		while(fedSize>=k){

			for(int i=0;i<fedSize-k;i++){
				//cerr<<"operating at "<<i<<endl;
				
				int absg0=curgstart0+i;
				curKmer=&buffer_c_str[i];
			
				/*if(StringUtil::toLower(buffer.substr(i,17))=="acaccataaatgcattc"){
					cerr<<"found at "<<(absg0+1)<<endl;	
				}*/
				
				if (strand==FORWARD_STRAND || strand==BOTH_STRANDS){
					//score forward strand here
					
					int score=scorer->score(curKmer);
					Hit thisHit(score,ref,absg0,absg0+k,FORWARD_STRAND);
					this->proposeHit(thisHit);
				}	
				if (strand==REVERSE_STRAND || strand==BOTH_STRANDS){
					//score reverse strand here
					
					int score=reverseScorer->score(curKmer);
					
					/*if (absg0==1919525){
						cerr<<"reverseScorer:"<<reverseScorer->toString()<<endl;
						cerr<<buffer.substr(i,17)<<" score at 1919526="<<score<<endl;	
					}*/
					
					Hit thisHit(score,ref,absg0,absg0+k,REVERSE_STRAND);
					this->proposeHit(thisHit);
				}	
				//cerr<<"done "<<i<<endl;
			}
			
			
			//advance curgstart0
			curgstart0+=fedSize-k+1;
			
			//now advance buffer window
			buffer=raf.get(curgstart0,min(curgstart0+bufferSize,gend1));
			buffer_c_str=buffer.c_str();
			fedSize=buffer.length();

		}
		
		
			
	}
	
	
};




int main(int argc,char** argv){
	/*Matrix<double> M;
	vector<double> row;
	row.push_back(0.0);
	row.push_back(0.5);
	row.push_back(0.7);
	row.push_back(0.8);
	
	for(int i=0;i<3;i++){
		for(int j=0;j<4;j++){
			row[j]+=1;	
		}	
		
		M.addRow(row);
	}
	
	M.print(cout);
	cerr<<"change 0,0"<<endl;
	M.itemAt(0,0)=-20;
	M.print(cout);
	M.finalize();
	cout<<"finalize"<<endl;
	M.print(cout);
	M.itemAt(1,2)=-10;
	cout<<"change 1,2"<<endl;
	M.print(cout);
	cout<<"read matrix"<<endl;
	Matrix<double> Mf;
	StringToFloat stf;
	Mf.readFile("matrix.dat",stf);
	Mf.print(cout);
	Mf.finalize();
	cout<<"finalize read matrix"<<endl;
	Mf.print(cout);*/
	
	/*RandomAccessFile hello("helloworld.txt");
	cout<<"1a:"<<hello.get(2,6)<<endl;
	cout<<"1b:"<<hello.get(4,9)<<endl;
	cout<<"1c:"<<hello.get(2,6)<<endl;
	cout<<"2:"<<hello.get(5,23)<<endl;
	//cout<<hello.get(100,12)<<endl;
	cout<<"3:"<<hello.get(2,6)<<endl;
	
	return 0;*/
	
	
	
	SimpleConsensusKmerScorer SryL("gaAtgcatttatggtgt",2);
	SimpleConsensusKmerScorer SryLR=SryL.createReverseComplementScorer();
	SimpleConsensusKmerScorer SryR("gctgggccaacttgtgcct",2);
	SimpleConsensusKmerScorer SryRR=SryR.createReverseComplementScorer();	
	/*string s;
	s="gaatgcatttatggtgt";cout<<s<<" "<<scks.score(s.c_str())<<endl;
	s="gaatgcacttatggtgt";cout<<s<<" "<<scks.score(s.c_str())<<endl;
	s="gaatgcacctatggtgt";cout<<s<<" "<<scks.score(s.c_str())<<endl;
	s="gaatgcacccatggtgt";cout<<s<<" "<<scks.score(s.c_str())<<endl;
	
	
	s=reverse_complement("gaatgcatttatggtgt");cout<<s<<" "<<scks.score(s.c_str())<<endl;
	s=reverse_complement("gaatgcatttatggtgt");cout<<s<<" "<<scksR.score(s.c_str())<<endl;
	s=reverse_complement("gaatgcacttatggtgt");cout<<s<<" "<<scksR.score(s.c_str())<<endl;
	s=reverse_complement("gaatgcacctatggtgt");cout<<s<<" "<<scksR.score(s.c_str())<<endl;
	s=reverse_complement("gaatgcacccatggtgt");cout<<s<<" "<<scksR.score(s.c_str())<<endl;	*/
	
	RandomAccessFile raf("/lab/jaenisch_albert/genomes/mm9/seq/chrY.seq");
	//s="gaatgcatttatggtgt";
	cout<<"SryL"<<endl;
	PWMHitsFinder SryLFinder(&SryL,&SryLR,20,PWMFINDER_BUFFER_SIZE);
	SryLFinder.scanSeq(raf,"chrY",0,numeric_limits<int>::max(),BOTH_STRANDS);
	SryLFinder.printHits(cout);
	cout<<"-------------------------------------------------------------------------"<<endl;
	
	//cout<<raf.get(1520,2530)<<endl;
	
	cout<<"SryR"<<endl;
	PWMHitsFinder SryRFinder(&SryR,&SryRR,20,PWMFINDER_BUFFER_SIZE);
	
	SryRFinder.scanSeq(raf,"chrY",0,numeric_limits<int>::max(),BOTH_STRANDS);
	SryRFinder.printHits(cout);
	
	raf.close();
	
}
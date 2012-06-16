#ifndef FASTAFILE_H_
#define FASTAFILE_H_

#include <fstream>
using namespace std;


class FastaFile
{
public:
	ifstream fin;
	string seqName;
	string seq;
	std::streampos last;
	bool upperCase;
	
	inline FastaFile(string filename,bool _upperCase=true): upperCase(_upperCase)
	{
		
		fin.open(filename.c_str());
		if(!fin.good())
		{
			fin.close();
			cerr<<"file "<<filename<<" cannot be opened"<<endl;
		}
	}
	
	inline bool readEntry()
	{
		seq="";
		seqName="";
		if(fin.eof())
			return false;
		
		fin>>seqName;

		if(seqName[0]!='>')
		{
			cerr<<"bad format!"<<endl;
			return false;
		}
		seqName=seqName.substr(1);
		last=fin.tellg();
		string line=" ";
		
		
		while(!fin.eof())
		{
			line="";
			fin>>line;
			//if(fin.eof())
			//	break;
			
			if(line[0]=='>')
			{
				fin.seekg(last);
				break;
				
			}
			
			if(upperCase)
				seq+=toUpperAndNoSpecialDNA(line);
			else
				seq+=line;
			last=fin.tellg();
		}
		
		return seq!="";
	}
	
	inline ~FastaFile()
	{
		fin.close();
	}
	
};

#endif /*FASTAFILE_H_*/

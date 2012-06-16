
#include <sam.h>
#include <string>
#include <set>
using namespace std;

#define BAMAUX_INTVALUE 'i'
#define BAMAUX_FLOATVALUE 'f'
#define BAMAUX_DOUBLEVALUE 'd'
#define BAMAUX_CHARVALUE 'A'
#define BAMAUX_STRINGVALUE 'Z'
#define BAMAUX_HAXSTRINGVALUE 'H'
#define BAMAUX_NOTSET '#'

class BamReader{
public:
	
	samfile_t*bamfile;
	bam_index_t*idx;
	string bamfilename;
	set<string> chromsInBam;
	
	class fetch_count_data{
	public:
		bam_header_t* bamHeader;
		int count;
		int start0Constraint;
		int end1Constraint;
		inline void resetCount(){
			count=0;
		}
		inline fetch_count_data(bam_header_t* _bamHeader,int _start0Constraint=-1,int _end1Constraint=-1,int _count=0):bamHeader(_bamHeader),count(_count),start0Constraint(_start0Constraint),end1Constraint(_end1Constraint){}
		
	};

	class fetch_count_data_divHits{
	public:
		bam_header_t* bamHeader;
		double count;
		int start0Constraint;
		int end1Constraint;
		inline void resetCount(){
			count=0.0;
		}
		inline fetch_count_data_divHits(bam_header_t* _bamHeader,int _start0Constraint=-1,int _end1Constraint=-1,double _count=0.0):bamHeader(_bamHeader),count(_count),start0Constraint(_start0Constraint),end1Constraint(_end1Constraint){}
		
	};	
	
	class BamAuxStruct{
	public:
		uint8_t type;
		union{
			uint8_t charValue;
			int64_t intValue;
			float floatValue;
			double doubleValue;
		};
		string stringValue;
		
		inline ostream& printToStream(ostream& os){
			switch(type){
				case BAMAUX_INTVALUE:
					os<<intValue;
					break;
				case BAMAUX_FLOATVALUE:
					os<<floatValue;
					break;
				case BAMAUX_DOUBLEVALUE:
					os<<doubleValue;
					break;
				case BAMAUX_CHARVALUE:
					os<<charValue;
					break;
				case BAMAUX_STRINGVALUE:
				case BAMAUX_HAXSTRINGVALUE:
					os<<stringValue;
					break;
				default:
					break;
			}
			
			return os;
		}
		
		inline BamAuxStruct():type(BAMAUX_NOTSET){
			//no need to initiliaze the number values, coz you should check type first. the stringValue is init to empty string by its default constructor anyway.
		}
		
	};
	
	inline void close(){
		if(bamfile){
			samclose(bamfile);	
			bamfile=NULL;	
		}
		if(idx){
			bam_index_destroy(idx);
			idx=NULL;
	
		}
	}
	
	inline ~BamReader(){
		close();
		
	}
	
	inline static int fetch_count(const bam1_t *b, void *data)
	{
		fetch_count_data* fdata=(fetch_count_data*)data;
		if(fdata->start0Constraint!=-1 && b->core.pos<fdata->start0Constraint)
			return 0;
		if(fdata->end1Constraint!=-1 && b->core.pos>=fdata->end1Constraint)
			return 0;
		fdata->count++;
		return 0;
	}
	
	inline static int fetch_count_divHits(const bam1_t *b, void *data)
	{
		fetch_count_data_divHits* fdata=(fetch_count_data_divHits*)data;
		if(fdata->start0Constraint!=-1 && b->core.pos<fdata->start0Constraint)
			return 0;
		if(fdata->end1Constraint!=-1 && b->core.pos>=fdata->end1Constraint)
			return 0;
		
		int numHits=BamReader::getNumHits(b,1); //if not given, assume num hits = 1
		
		
		fdata->count+=1.0/numHits;
		return 0;
	}	
	
	inline BamReader(const string& _bamfilename):bamfilename(_bamfilename),idx(NULL),bamfile(NULL){
		
		bamfile=samopen(bamfilename.c_str(),"rb",0);
		
		if(!bamfile)
		{
			cerr<<"Fail to open BAM file "<<bamfilename<<endl;
			//return 1;
		}
		
		cerr<<"load bamfile index"<<endl;
		
		idx=bam_index_load(bamfilename.c_str());
		if(!idx){
			cerr<<"Fail to load index\n"<<endl;
			//return 1;
		}
		
		for(int i=0;i<bamfile->header->n_targets;i++){
			cerr<<"bam file has reference "<<bamfile->header->target_name[i]<<endl;
			chromsInBam.insert(bamfile->header->target_name[i]);
		}
		
	}
	
	inline static string getQName(bam1_t* bamT){
		char qnameBuffer[256];
		strncpy(qnameBuffer,bam1_qname(bamT),bamT->core.l_qname);
		return qnameBuffer;	
		
	}
	
	inline static bool hasMultipleFragments(bam1_t* bamT){
		return (bamT->core.flag & 0x1);
	}
	
	inline static bool isFirstFragment(bam1_t* bamT){
		return (bamT->core.flag & 0x40);
	}
	
	inline static bool eachFragmentProperlyAligned(bam1_t* bamT){
		return (bamT->core.flag & 0x2);
	}
	
	inline static bool fragmentUnmapped(bam1_t* bamT){
		return (bamT->core.flag & 0x4);
	}
	
	inline static void getAuxFields(const bam1_t*b,map<string,BamAuxStruct>&  auxfields,bool clearMap=true){
		uint8_t* s=bam1_aux(b);
		
		if(clearMap)
			auxfields.clear();
		
		//copied and modified from bam.c from samtools
		while (s < b->data + b->data_len) {
			uint8_t type, key[3];
			key[0] = s[0]; key[1] = s[1]; key[2]='\0'; //add in terminator char '\0'
			s += 2; type = *s; ++s;
			
			string skey((char*)key);
			BamAuxStruct bas;
			
			
			
			if (type == 'A') { bas.type=BAMAUX_CHARVALUE; bas.charValue=*s; ++s; }
			else if (type == 'C') { bas.type=BAMAUX_INTVALUE; bas.intValue=int64_t(*s); ++s; }
			else if (type == 'c') { bas.type=BAMAUX_INTVALUE; bas.intValue=int64_t(*(int8_t*)s); ++s; }
			else if (type == 'S') { bas.type=BAMAUX_INTVALUE; bas.intValue=int64_t(*(uint16_t*)s); s += 2; }
			else if (type == 's') { bas.type=BAMAUX_INTVALUE; bas.intValue=int64_t(*(int16_t*)s); s += 2; }
			else if (type == 'I') { bas.type=BAMAUX_INTVALUE; bas.intValue=int64_t(*(uint32_t*)s); s += 4; }
			else if (type == 'i') { bas.type=BAMAUX_INTVALUE; bas.intValue=int64_t(*(int32_t*)s); s += 4; }
			else if (type == 'f') { bas.type=BAMAUX_FLOATVALUE;bas.floatValue=(*(float*)s); s += 4; }
			else if (type == 'd') { bas.type=BAMAUX_DOUBLEVALUE; bas.floatValue=(*(double*)s); s += 8; }
			else if (type == 'Z' || type == 'H') { bas.type=type; while (*s) bas.stringValue+=(*s++); ++s; }
	
			auxfields.insert(map<string,BamAuxStruct>::value_type(skey,bas));
		}		
		
	}
	
	inline static BamAuxStruct getAuxFieldValue(const bam1_t*b,const string& targetKey){
		uint8_t* s=bam1_aux(b);
		
		BamAuxStruct bas;
		//copied and modified from bam.c from samtools
		while (s < b->data + b->data_len) {
			uint8_t type, key[3];
			key[0] = s[0]; key[1] = s[1]; key[2]='\0'; //add in terminator char '\0'
			s += 2; type = *s; ++s;
			
			string skey((char*)key);
			
			
			bool matched=(skey==targetKey);
			
			if (type == 'A') { if(matched){bas.type=BAMAUX_CHARVALUE; bas.charValue=*s; return bas;} ++s; }
			else if (type == 'C') { if(matched){bas.type=BAMAUX_INTVALUE; bas.intValue=int64_t(*s); return bas;} ++s; }
			else if (type == 'c') { if(matched){bas.type=BAMAUX_INTVALUE; bas.intValue=int64_t(*(int8_t*)s); return bas;} ++s; }
			else if (type == 'S') { if(matched){bas.type=BAMAUX_INTVALUE; bas.intValue=int64_t(*(uint16_t*)s); return bas;} s += 2; }
			else if (type == 's') { if(matched){bas.type=BAMAUX_INTVALUE; bas.intValue=int64_t(*(int16_t*)s); return bas;} s += 2; }
			else if (type == 'I') { if(matched){bas.type=BAMAUX_INTVALUE; bas.intValue=int64_t(*(uint32_t*)s); return bas;} s += 4; }
			else if (type == 'i') { if(matched){bas.type=BAMAUX_INTVALUE; bas.intValue=int64_t(*(int32_t*)s); return bas;} s += 4; }
			else if (type == 'f') { if(matched){ bas.type=BAMAUX_FLOATVALUE;bas.floatValue=(*(float*)s); return bas;}s += 4; }
			else if (type == 'd') { if(matched){bas.type=BAMAUX_DOUBLEVALUE; bas.floatValue=(*(double*)s); return bas;} s += 8; }
			else if (type == 'Z' || type == 'H') { if(matched){bas.type=type; while (*s) bas.stringValue+=(*s++);} else{ while(*s) s++;} ++s; }
			
			
		}	
		
		return bas; /* can have type BAMAUX_NOTSET */
		
	}
	
	inline static bool hasAuxField(const bam1_t* b,const string& targetKey){
		BamAuxStruct bas=getAuxFieldValue(b,targetKey);
		return (bas.type!=BAMAUX_NOTSET);
	}
	
	inline static void appendAuxField(bam1_t*b,const string& stag, BamAuxStruct* bas){
		
		const char* cstag=stag.c_str();
		
		char tag[2];
		tag[0]=cstag[0];
		tag[1]=cstag[1];
		
		switch(bas->type){
			case BAMAUX_INTVALUE:
				//assume signed int32
				bam_aux_append(b,tag,bas->type,sizeof(int32_t),(uint8_t*)&bas->intValue);
				break;
			case BAMAUX_FLOATVALUE:
				bam_aux_append(b,tag,bas->type,sizeof(float),(uint8_t*)&bas->floatValue);
				break;
			case BAMAUX_DOUBLEVALUE:
				bam_aux_append(b,tag,bas->type,sizeof(float),(uint8_t*)&bas->doubleValue);
				break;
			case BAMAUX_STRINGVALUE:
			case BAMAUX_HAXSTRINGVALUE:
				bam_aux_append(b,tag,bas->type,bas->stringValue.length(),(uint8_t*)bas->stringValue.c_str());
				break;
			case BAMAUX_CHARVALUE:
				bam_aux_append(b,tag,bas->type,sizeof(uint8_t),(uint8_t*)&bas->charValue);
			default:
				//do nothing don't know what to do
				break;
		}
	}
	
	inline static uint32_t getNumHits(const bam1_t*bamT,uint32_t defaultValue){
		BamAuxStruct bas=getAuxFieldValue(bamT,"NH");
		if(bas.type==BAMAUX_INTVALUE){
			return bas.intValue;
		}
			
		return defaultValue;
		
		
	}
	

	
	inline static bool fragmentMapped(bam1_t* bamT){
		return !fragmentUnmapped(bamT);
	}
	
	inline static bool nextFragmentUnmapped(bam1_t*bamT){
		return (bamT->core.flag & 0x8);
	}
	
	inline static bool nextFragmentMapped(bam1_t* bamT){
		return !nextFragmentUnmapped(bamT);
	}
	inline static bool seqReverseComplemented(bam1_t*bamT){
		return (bamT->core.flag & 0x10);
	}
	inline static bool seqNextFragmentReverse(bam1_t*bamT){
		return (bamT->core.flag & 0x20);
	}
	
	inline static bool isLastFragment(bam1_t* bamT){
		return (bamT->core.flag & 0x80);
	}
	
	inline static bool secondaryAlignment(bam1_t*bamT){
		return (bamT->core.flag & 0x100);
	}
	
	inline static bool notPassQC(bam1_t*bamT){
		return (bamT->core.flag & 0x200);
	}
	
	inline static bool PCROrOpticalDuplicate(bam1_t*bamT){
		return (bamT->core.flag & 0x400);
	}
	
	inline static unsigned int naiveFragmentCount(bam1_t*bamT){
		if(fragmentMapped(bamT)){
			if(hasMultipleFragments(bamT))
			{
				if(isFirstFragment(bamT))
				{
					return 1;
				}
				else 
				{
					if(nextFragmentUnmapped(bamT))
					{
						return 1;
					}
				}

			}
			else
			{
				return 1;
			}
		}
		
		return 0;
	}
	
	class FragmentSetCounter{
	public:
		set<string> qnames;
		inline void countRead(bam1_t*bamT){
			qnames.insert(getQName(bamT));
		}
		inline operator unsigned int(){
			return qnames.size();
		}
	};
	
	inline static string getGenomeBrowserCoordinate(const string&chrom,int start0,int end1){
		char buf[256];
		sprintf(buf,"%s:%d-%d",chrom.c_str(),start0+1,end1);
		return buf;
	}
	
	inline static unsigned char getMappingQualityScore(bam1_t*bamT){
		return bamT->core.qual;	
	}
	
	inline static int countTotalNumOfReadsInBam(string bamfile){
		int total=0;
		samfile_t* bf=samopen(bamfile.c_str(),"rb",0);
		
		
		
		if(!bf){
			cerr<<"bam file "<<bamfile<<" cannot be open for counting"<<endl;
			return -1;
		}
		
		//now do the counting
		bam1_t *dummy=bam_init1();
		
		
		while(samread(bf,dummy)>=0){
			total++;
			if(total%1000000==1){
				cerr<<"passing through read "<<total<<endl;
			}
		}
		
		bam_destroy1(dummy);
		samclose(bf);
		return total;
	}
	
	inline int fetchCountOverlappingRegion(const string&chrom, int start0, int end1,bool constrainStart){
		string region=getGenomeBrowserCoordinate(chrom,start0,end1);
		int ref_id;
		int begin;
		int end;
		
		fetch_count_data counter(bamfile->header,(constrainStart?start0:-1));
		
		bam_parse_region(bamfile->header,region.c_str(),&ref_id,&begin,&end);
		bam_fetch(bamfile->x.bam,idx,ref_id,begin,end,&counter,BamReader::fetch_count);
		
		return counter.count;
	}
	
	inline double fetchCountOverlappingRegionDivideByNumHits(const string&chrom, int start0, int end1,bool constrainStart){
		string region=getGenomeBrowserCoordinate(chrom,start0,end1);
		int ref_id;
		int begin;
		int end;
		
		fetch_count_data_divHits counter(bamfile->header,(constrainStart?start0:-1));
		
		bam_parse_region(bamfile->header,region.c_str(),&ref_id,&begin,&end);
		bam_fetch(bamfile->x.bam,idx,ref_id,begin,end,&counter,BamReader::fetch_count_divHits);
		
		return counter.count;
	}
	
};

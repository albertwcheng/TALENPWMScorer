/***************************************************************************
 Copyright 2010 Wu Albert Cheng <albertwcheng@gmail.com>
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
*******************************************************************************/ 

#ifndef NUCLEIC_H_
#define NUCLEIC_H_

#define A_SYM			'A'
#define A_SYML			'a'
#define C_SYM			'C'
#define C_SYML			'c'
#define G_SYM			'G'	
#define G_SYML			'g'
#define T_SYM			'T'
#define T_SYML			't'	
#define U_SYM			'U'
#define U_SYML			'u'
#define N_SYM			'N'
#define N_SYML			'n'
#define X_SYM			'X'
#define X_SYML			'x'
#define R_SYM			'R'
#define R_SYML			'r'
#define Y_SYM			'Y'
#define Y_SYML			'y'
#define END_SYM			'.'
#define END_SYML		END_SYM




#define GAP_SYM1		'-'
#define GAP_SYM2		'.'
#define GAP_SYM3		' '


#define BRACKET			'('
#define CLOSE_BRACKET 	')'

#define A_IND	0
#define C_IND	1
#define G_IND	2
#define T_IND	3
#define N_IND	4
#define X_IND	5

#include <string.h>


inline int getIndex(char c)
{
	switch(c)
	{
	case A_SYM:case A_SYML:
		return A_IND;
	case T_SYM:case T_SYML:
		return T_IND;
	case C_SYM:case C_SYML:
		return C_IND;
	case G_SYM:case G_SYML:
		return G_IND;
	case U_SYM:case U_SYML:
		return T_IND;
	default:
		return N_IND;
	}	
}

inline char complementBase(char c)
{
	switch(c)
	{
	case A_SYM:case A_SYML:
		return T_SYM;
	case T_SYM:case T_SYML:
		return A_SYM;
	case C_SYM:case C_SYML:
		return G_SYM;
	case G_SYM:case G_SYML:
		return C_SYM;
	case U_SYM:case U_SYML:
		return A_SYM;
	default:
		return c;
	}	
}
inline string complement(const string& src)
{
	string dst;
	
	const char* srcs=src.c_str();
	int len=src.length();
	char* dsts=new char[len+1];
	
	int j=0;
	for(int i=0;i<len;i++,j++)
	{
		dsts[j]=complementBase(srcs[i]);
	}
	
	dsts[j]='\0';
	dst.append(dsts);
	delete[] dsts;
	return dst;
}

inline string reverse(const string& src)
{
	string dst;
	
	const char* srcs=src.c_str();
	int len=src.length();
	char* dsts=new char[len+1];
	
	int j=0;
	for(int i=len-1;i>=0;i--,j++)
	{
		dsts[j]=srcs[i];
	}
	
	dsts[j]='\0';
	dst.append(dsts);
	delete[] dsts;
	return dst;
}

inline string reverse_complement(const string& src)
{
	string dst;
	
	const char* srcs=src.c_str();
	int len=src.length();
	char* dsts=new char[len+1];
	
	int j=0;
	for(int i=len-1;i>=0;i--,j++)
	{
		dsts[j]=complementBase(srcs[i]);
	}
	
	dsts[j]='\0';
	dst.append(dsts);
	delete[] dsts;
	return dst;
}

#endif /*NUCLEIC_H_*/

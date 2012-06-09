#ifndef _SYSTEM_UTIL_H_
#define _SYSTEM_UTIL_H_
/***************************************************************************
 Copyright 2011 Wu Albert Cheng <albertwcheng@gmail.com>
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
 
#include <libgen.h>
#include <string.h>
#include <StringUtil.h>

#define PATH_SEPARATOR "/"

class SystemUtil
{
	public:
		inline static string basename(const string& path)
		{
			return basename(path.c_str());
			
		}	
		inline static string dirname(const string& path)
		{
			char*tmpString=strdup(path.c_str());
			string returna(dirname(tmpString));
			delete[] tmpString;
			return returna;
			
			
		}
		
		inline static string formPathFromComponents(const vector<string>& components){
			return StringUtil::join<vector<string>,vector<string>::const_iterator>(components,PATH_SEPARATOR);
		}
};

#endif /*_SYSTEM_UTIL_H*/


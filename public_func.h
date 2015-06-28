
#ifndef _CORE_PUBLIC_FUNC_H_
#define _CORE_PUBLIC_FUNC_H_

#include<string>
#include<sstream>

class PubFuncs
{
public:
	static std::string cvtInt2Str(int i); //Convert int to string 
		
	static int cvtStr2Int(std::string str); //Convert string to int 

	//static bool fileExist(std::string fpath);// Check whether a file exist or not 
};

#endif
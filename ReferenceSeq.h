#ifndef _H_ECNEREFER_
#define _H_ECNEREFER_

#include<string>

#include"fasta_parser.h"

class ReferenceSeq
{
public:
	ReferenceSeq();
	ReferenceSeq(std::string path);

public:
	int getLength();//return the length of the reference.
	void getSubRef(std::string chrom,int pos, int length, std::string& sub_ref);//return sub-reference. 
	void getSupplementSubRef(std::string chrom,int pos, int length, std::string& ssuplmt);//return supplement sub-reference. 
	void getSupplementSubRef(std::string& sub_ref, std::string& ssuplmt);

	void mapChromIDName();//map chrom ID with Name, by using XX.fasta.fai file 
	void loadChromIDName();
	void getChromIDByName(std::string name, int& id);//given chrom name, return chrom id.
	void getChromNameByID(int id, std::string& name);//given chrom id, return chrom name.

public:
	void setPath(std::string path);

private:
	void cvt2Uppercase();//convert Caps to Lower case 
	void supplementRef(std::string& sorgnl, std::string& ssuplmt);//convert to summplement ref 

private:
	std::string path;//path of the reference sequence 
	FastaParser fp;
};

#endif

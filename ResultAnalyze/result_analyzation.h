/************************************************************
Description:
	Calc the accurary of called exons.
    true_exons and called_exons both should be sorted by chrom and pos
*************************************************************/


#ifndef _H_NOITAZYLANA_TLUSER_
#define _H_NOITAZYLANA_TLUSER_

#include<string>
#include<vector>
#include"./../exon.h"

class RsltAnalyzation
{
public:
	RsltAnalyzation();
	RsltAnalyzation(std::string texons, std::string cexons, int offset);

public:
	void cmpExonAccuracy(double& accuracy);//calc the accuracy by comapring with true exons 

public:
	void setFlex(int iflex);

private:
	void loadExons(std::string filename, std::vector<Exon>& exons);//load exons from file into memory

private:
	std::string sfile_true_exons;//path of the true exons file.
	std::string sfile_called_exons;//path of the called exons file.
	int iflex;// flexible value
};

#endif

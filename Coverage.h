/*
Coverage.h
Description:
	Calc the reads coverage, and classify the reads. 
*/

#ifndef _H_EGAREVOC_
#define _H_EGAREVOC_

#include<string>
#include<vector>
#include"bam_parse.h"
#include"CandidateSitesCaller.h"

class Coverage
{
public:
	Coverage();
	Coverage(std::string bam_name, std::string fai_path);
	std::string getBamPath();
	void setBamPath(std::string bam_name);
	void setInsertSize(double mean, double devi);
	void initFiles();

public:
	void pileup(int cid, BamParse& bp);

private:
	void pileupRegion(int cid, int start, int end, BamParse& bp, int*& bases, Brkpnt*& brkpnts, int*& is_clip_pos, int*& ldiscor, int*& rdiscor,\
	ofstream& fout_fm, ofstream& fout_clip, ofstream& fout_um,ofstream& fout_others, ofstream& fout_cov,\
	ofstream& fout_lclip, ofstream& fout_rclip, string& chrom_name);
	
	void callCandidateByCov(Brkpnt*& brkpnts, int*& cov, int*& is_clip_pos, int*& ldiscor, int*& rdiscor, int chrom_len);
	int isDiscordant(int rid, int map_pos, int rnext, int pnext, int discor_range, int max_range);
	void saveBrkpntInfo(int chrom_len, int*& is_clip_pos, int*& ldiscor, int*& rdiscor);

private:
	std::string bam_path;

private:
	double mean_insert;//mean insert size
	double devi_insert;//standard variaiton of insert size
};

#endif


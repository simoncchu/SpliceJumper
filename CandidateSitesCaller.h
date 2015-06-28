#ifndef _H_RELLACSETISETADIDNAC_
#define _H_RELLACSETISETADIDNAC_

#include<string>
#include<vector>
#include<utility>
#include"bam_parse.h"
#include"ReferenceSeq.h"

using namespace std;

class Brkpnt
{
public:
	Brkpnt(){pos=0; lclip=0; rclip=0;}
	Brkpnt(int pos, int lclip, int rclip){this->pos=pos; this->lclip=lclip; this->rclip=rclip;}
	int pos;
	int lclip;
	int rclip;
};


class JunctionFeatures
{
public:
	JunctionFeatures(){}
	JunctionFeatures(int start, int end, int softclip_map, int iscanon, double cov_diff, int full_map, int discor) 
	{this->start=start; this->end=end; this->softclip_map=softclip_map, this->iscanon=iscanon; \
	this->cov_diff=cov_diff; this->full_map=full_map; this->discor=discor;}

public:
	bool operator < (const JunctionFeatures& jf) const
	{
		if(this->start==jf.start)
			return this->end < jf.end;
		else
			return this->start <  jf.start;
	}

public:
	int start;
	int end;
	int softclip_map;
	int iscanon;
	double cov_diff;
	int full_map;
	int discor;
};

class BrkConnectNode
{
public:
	int id;
	int pos;
	int split_support;
	int pe_support;
	BrkConnectNode* pnext;
};

class SiteFeatures 
{
public:
	SiteFeatures(){softclip_map=0; lclip=0; rclip=0; lcov=0.0; rcov=0.0; bdonor=true; vconnct.clear();}
	SiteFeatures(int pos, int lclip, int rclip, bool bdonor){this->pos=pos; this->lclip=lclip; this->rclip=rclip;\
		softclip_map=0; lcov=0.0; rcov=0.0; this->bdonor=bdonor; vconnct.clear();}

public:
	int pos;//site position
	int softclip_map; //soft-clipped segment can be mapped here
	int lclip; //left clipped reads
	int rclip;//right clipped reads
	//int full_map;//# of fully mapped reads
	double lcov; //left region coverage
	double rcov;//right region coverage
	bool bdonor;//is dornor site
	int discor;//# of discordant pairs encompassing the brkpnt
	std::vector<std::pair<int,int> > vconnct; // < connect site index in vbrkpnts[], # of supported reads> 
	std::vector<bool> vcanon;//whether is canonical signal for each pair <site, vconnct[i]>
};


class CandidateSitesCaller
{
public:
	CandidateSitesCaller();
	~CandidateSitesCaller();
	void setRefPath(std::string fref);
	void setChromNameLen(std::string chrom_name, int chrom_len);
	void setInsertSize(double mean, double devi);

public:
	void callCandidateSites();
	void outputFeatures(std::string fout_path);

	void outputGraph();
	void releaseGraph();

private:
	void loadBrkpntInfo(int*& bases, int*& is_clip_pos, int*& ldiscor, int*& rdiscor);//load information from files
	void alignSoftClip(std::string frefs, std::string freads, int*& is_clip_pos, int*& clip_pos_order);
	void getCigarPosLen(std::vector< std::pair<std::string,int> >& cigar, int seg_len, int& pos_front, int& len_front, int& pos_tail, int& len_tail);
	bool searchMapSite(std::string segmnt, ReferenceSeq& refseq, std::string chrom_name, int vi, int pnext, bool bfront);
	bool alignSegmt(std::string sub_ref, std::string str_clip);
	bool isDonorSite(string str, std::size_t& offset);
	bool isAcceptorSite(string str, std::size_t& offset);
	void loadHardClipRawReads(map<string, string>& lhclip_map, map<string, string>& rhclip_map);

public:
	std::vector<SiteFeatures> vbrkpnts;
	BrkConnectNode** brkcon;
	int nbrkpnt;

private:
	std::string fref;// full path of the reference file
	std::string chrom_name;//chrom name 
	int chrom_len;//chrom length

	double mean_insert;//mean insert size
	double devi_insert;//standard variaiton of insert size
};

#endif 

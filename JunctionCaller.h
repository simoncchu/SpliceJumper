#ifndef _H_RELLACNOITCNUJ_
#define _H_RELLACNOITCNUJ_

#include<vector>
#include<utility>
#include<string>
#include"CandidateSitesCaller.h"

class JunctionCaller
{
public:
	JunctionCaller(){vbrkpnts.clear();}
	~JunctionCaller();
	void setFiles(std::string fbrkpnt, std::string fjunctions);
	
public:
	void setChromLen(int chrom_len);

public:
	void prepareBrkpnts();
	void callJunctions(); 
	
private:
	bool isTrueSite(std::vector<pair<int,bool> >& vbrkpnts, int pos);
	void callJunctionsByPEReads();
	//void loadGraph(std::string fgraph); //load brkpnt relationship graph

private:
	std::string fbrkpnt; //file save called out sites
	std::string fjunctions;//file save all the candidate junctions
	std::vector<std::pair<int,bool> > vbrkpnts;
	int chrom_len;
	std::vector<std::pair<int,int> > vconnct;
	BrkConnectNode** brkcon;
};

#endif

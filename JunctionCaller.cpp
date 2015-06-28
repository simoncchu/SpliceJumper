#include"JunctionCaller.h"
#include<fstream>

using namespace std;

JunctionCaller::~JunctionCaller()
{
	vector<pair<int,bool> >().swap(this->vbrkpnts);
}

void JunctionCaller::setFiles(std::string fbrkpnt, std::string fjunctions)
{
	this->fbrkpnt=fbrkpnt;
	this->fjunctions=fjunctions;
}

void JunctionCaller::prepareBrkpnts()
{
	ifstream fin_brkpnt;
	fin_brkpnt.open(fbrkpnt.c_str());//file save all the true brkpnts
	vbrkpnts.clear();
	int pos;
	int idonor;
	bool bdonor;
	while(fin_brkpnt>>pos>>idonor)
	{
		if(idonor==1) bdonor=true;
		else bdonor=false;
		vbrkpnts.push_back(std::make_pair(pos,bdonor));
	}
	fin_brkpnt.close();
}


void JunctionCaller::callJunctions()
{
	//last_start<<" "<<last_end<<" "<<last_support<<" "<<last_fullmap<<" "<<last_discor<<" "<<last_covdiff<<" "<<last_iscanon
	ofstream fout_jct; //output called out junction file 
	fout_jct.open("junctions.txt");

	int pos1, pos2;
	ifstream fin_cj;//file save all the candidate junctions ifstream
	fin_cj.open(this->fjunctions.c_str());
	int start,end, support,fullmap,discor,iscanon;
	double covdiff;
	while(fin_cj>>start>>end>>support>>fullmap>>discor>>covdiff>>iscanon)
	{
		//check whether start and end are exist in vbrkpnts
		if((isTrueSite(vbrkpnts, start)==true) && (isTrueSite(vbrkpnts,end)) )
		{
			fout_jct<<start<<" "<<end<<endl;
		}
	}

	fout_jct.close();
	fin_cj.close();
	
}


//----------private functions---------------------------------------------------------------------------------
void JunctionCaller::callJunctionsByPEReads()
{
	int* plus1=new int[this->chrom_len];
	memset(plus1,-1,chrom_len);
	int* minus1=new int[this->chrom_len];
	memset(minus1,-1,chrom_len);

	int nbrkpnt=this->vbrkpnts.size();
	for(int i=0;i<nbrkpnt;i++)
	{
		plus1[vbrkpnts[i].first]=0;
		minus1[vbrkpnts[i].first]=0;
	}

	int pcnt=-1;
	int mcnt=0;
	for(int i=0; i<this->chrom_len; i++)
	{
		if(minus1[i]==0)
		{
			mcnt++;
			minus1[i]=mcnt;
		}
		else
		{
			minus1[i]=mcnt;
		}

		if(plus1[i]==0)
		{
			plus1[i]=pcnt;
			pcnt++;
		}
		else
		{
			plus1[i]=pcnt;
		}
	}

	//Generate another graph
	//check each reads
	//first check each fully mapped reads

	//then check each clipped mapped reads
	
	//

	delete[] plus1;
	delete[] minus1;
}

bool JunctionCaller::isTrueSite(vector<pair<int,bool> >& vbrkpnts, int pos)
{
	int vsize=vbrkpnts.size();
	bool bexist=false;
	for(int i=0;i<vsize;i++)
	{
		if(pos==vbrkpnts[i].first)
		{
			bexist=true;
			break;
		}
	}
	return bexist;
}

//------------------------------------------------------------------------------------------------------
void  JunctionCaller::setChromLen(int chrom_len)
{
	this->chrom_len=chrom_len;
}
#include"result_analyzation.h"
#include<fstream>
#include <cmath>
#include<iostream>

using namespace std;

RsltAnalyzation::RsltAnalyzation()
{
	this->iflex=50;
}

RsltAnalyzation::RsltAnalyzation(std::string texons, std::string cexons, int iflex)
{
	this->sfile_true_exons=texons;
	this->sfile_called_exons=cexons;
	this->iflex=iflex;
}

void RsltAnalyzation::cmpExonAccuracy(double& accuracy)
{
	std::vector<Exon> texons;//true exons
	std::vector<Exon> cexons;//called exons 
	texons.clear();
	cexons.clear();

	loadExons(this->sfile_true_exons, texons);
	loadExons(this->sfile_called_exons, cexons);

	int tsize=texons.size();
	int csize=cexons.size();
	
	cout<<"No. of true exons are:"<<tsize<<endl;
	cout<<"No. of called exons are:"<<csize<<endl;

	int t=0;//pointer to texons[t]
	int c=0;//pointer to cexons[c]

	int tchrom,tstart,tend;
	int cchrom,cstart,cend;

	int cnt_hit=0;
	while(t<tsize && c<csize)
	{
		tchrom=texons[t].getChrom();
		cchrom=cexons[c].getChrom();
		tstart=texons[t].getStartPos();
		cstart=cexons[c].getStartPos();
		tend=texons[t].getEndPos();
		cend=cexons[c].getEndPos();

		if((tchrom > cchrom) || (tstart-cstart>iflex))
		{
			c++;
		}
		else if((tchrom < cchrom) || (cstart-tstart)>iflex) 
		{
			t++;
		}
		else
		{//same chromosome && start at the acceptable range

			if(abs(tstart-cstart)<=iflex && abs(tend-cend)<=iflex)
			{
				cnt_hit++;
			}

			if(tstart>cstart)
				c++;
			else if(cstart>tstart)
				t++;
			else
			{
				c++;
				t++;
			}

		}//end of else
	}//end of while

	cout<<"No. of hit exons is: "<<cnt_hit<<endl;

	accuracy=(double)cnt_hit/(double)tsize;
	
}

void RsltAnalyzation::setFlex(int iflex)
{
	this->iflex=iflex;
}

//----------------------------private functions-------------------------------------------------------------------------

//load exons from file into memory
void RsltAnalyzation::loadExons(std::string filename, std::vector<Exon>& exons)//load exons from file into memory
{
	ifstream fexon;
	fexon.open(filename.c_str(), ifstream::in);
	int chrom, start,end;
	while(fexon>>chrom>>start>>end)
	{
		Exon exon(chrom,start,end);
		exons.push_back(exon);
	}
	fexon.close();
}



/*

int exonStatistic()
{
	ifstream fin;
	fin.open("candiate_exon_all.site",ifstream::in);

	ofstream fout;
	fout.open("unique_candidate_exon.site", ofstream::out);
	
	int chromID=-1, pos=-1, direction=-1;
	fin>>chromID>>pos>>direction;
	int tempID=chromID,tempPos=pos,tempDirection=direction;
	int cnt_f=0, cnt_b=0;
	while(fin>>chromID>>pos>>direction)
	{
		if(tempID != chromID || tempPos != pos)
		{
			fout<<tempID<<" "<<tempPos<<" "<<cnt_f<<" "<<cnt_b<<endl;
			if(direction==1)
			{
				cnt_f=1;
				cnt_b=0;
			}
			else
			{
				cnt_b=1;
				cnt_f=0;
			}
			tempID=chromID;
			tempPos=pos;
		}
		else if(tempID==chromID && tempPos==pos)
		{
			if(direction==1)
				cnt_f+=1;
			else
				cnt_b+=1;
		}
		
	}

	fout.close();
	fin.close();
}

void exonAnalyze()
{
	ifstream fin;
	fin.open("unique_candidate_exon.site", ifstream::in);

	//ofstream fout;
	//fout.open("called_candidate_exon.site", ifstream::out);

	int chromID, pos, cntf,cntb;

	int cnt_ff=0;
	int cnt_bb=0;
	int cnt_exon=0;
	int cnt_all_sites=0;

	bool bpositive=false;
	bool bpre=false;
	int prechrom=-1;
	int prepos=-1;
	while(fin>>chromID>>pos>>cntf>>cntb)
	{
		cnt_all_sites++;
		
		int cov=cntf+cntb;
		if(cov<=5) continue;
		if(cntf>cntb)
			bpositive=true;
		else
			bpositive=false;
		
		if(chromID != prechrom)
		{
			continue;
		}

		if(bpre==true && bpositive==false)
		{
			//fout<<
			cnt_exon++;
		}
		else if(bpre==true && bpositive==true)
		{
			cnt_ff++;
		}
		else if(bpre==false && bpositive==false)
		{
			cnt_bb++;
		}
		bpre=bpositive;
	}
	cout<<" "<<cnt_all_sites<<" "<<cnt_exon<<" "<<cnt_ff<<" "<<cnt_bb<<endl;

	//fout.close();
	fin.close();
}

*/
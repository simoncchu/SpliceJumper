#include<string>
#include<vector>
#include<fstream>
#include<vector>
#include<utility>
#include<cmath>
#include<iostream>
#include<algorithm>
#include<map>

#include"CandidateSitesCaller.h"
#include"bam_parse.h"
#include"Coverage.h"
#include"public_parameters.h"
#include"ReferenceSeq.h"
#include"Alignment.h"
#include"./algorithms/local_alignment.h"
using namespace std;

CandidateSitesCaller::CandidateSitesCaller()
{
	vbrkpnts.clear();
}

CandidateSitesCaller::~CandidateSitesCaller()
{
	vector<SiteFeatures>().swap(this->vbrkpnts);
	this->releaseGraph();
}

void CandidateSitesCaller::setChromNameLen(std::string chrom_name, int chrom_len)
{
	this->chrom_len=chrom_len;
	this->chrom_name=chrom_name;
}


//load brkpnt information, all these information are saved when calc coverage
void CandidateSitesCaller::loadBrkpntInfo(int*& bases, int*& is_clip_pos, int*& ldiscor, int*& rdiscor)
{
	//load in coverage
	ifstream  fin_cov;
	fin_cov.open(fname_cov.c_str());
	string schrom;
	int index, icov;
	while(fin_cov>>schrom>>index>>icov)
	{
		bases[index]=icov;
	}
	fin_cov.close();

	//load in candidate brkpnts
	ifstream fin_cs;
	fin_cs.open(fname_brkpnts.c_str());
	string sdir;
	int pos, lclip, rclip;
	int clip;
	while(fin_cs>>chrom_name>>pos>>lclip>>rclip>>sdir)//read in candidate sites 
	{
		bool bdonor=false;
		SiteFeatures temp_brkpnt(pos, lclip, rclip, bdonor);
		
		if(sdir=="+") temp_brkpnt.bdonor=false;
		else if(sdir=="-") temp_brkpnt.bdonor=true;
		else
		{//sdir=="*"
			temp_brkpnt.lcov=-1.0;
			temp_brkpnt.rcov=-1.0;
		}
		vbrkpnts.push_back(temp_brkpnt);
	}
	fin_cs.close();

	//load in whether is clip pos
	ifstream fin_icp;
	fin_icp.open(fname_icp.c_str());
	int icp;
	while(fin_icp>>index>>icp)
	{
		is_clip_pos[index]=icp;
	}
	fin_icp.close();

	//load in ldiscor
	ifstream fin_ld;
	fin_ld.open(fname_ldiscor.c_str());
	int ileft;
	while(fin_ld>>index>>ileft)
	{
		ldiscor[index]=ileft;
	}
	fin_ld.close();

	//load in rdiscor
	ifstream fin_rd;
	fin_rd.open(fname_rdiscor.c_str());
	int iright;
	while(fin_rd>>index>>iright)
	{
		rdiscor[index]=iright;
	}
	fin_rd.close();
}

//void CandidateSitesCaller::callCandidateSites(BamParse& bp, int cid)
void CandidateSitesCaller::callCandidateSites()
{
	int* bases=new int[chrom_len+1];//save the coverage of each chromosome
	int* is_clip_pos=new int[chrom_len+1]; //record whether this position is an candidate position, combined ones record the distance to the record one.
	int* ldiscor=new int[chrom_len+1];
	int* rdiscor=new int[chrom_len+1];
	for(int i=0;i<chrom_len;i++)
	{
		bases[i]=0;
		is_clip_pos[i]=NO_CLIP_POS;
		ldiscor[i]=0;
		rdiscor[i]=0;
	}

	loadBrkpntInfo(bases, is_clip_pos, ldiscor, rdiscor);

	cout<<"Load brkpnt information successfully--------------------------------------------------------------"<<endl;//////////////////////////////////////////////////
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//ofstream fout_isclip;//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//fout_isclip.open("test_is_clip.txt"); 
	//for(int i=0;i<chrom_len;i++)
	//{
	//	if(is_clip_pos[i]!=NO_CLIP_POS)
	//	{
	//		fout_isclip<<i<<" "<<is_clip_pos[i]<<endl;
	//	}
	//}
	//fout_isclip.close();//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	string sdir;
	int pos, lclip, rclip;
	int clip;
	int discor_range=mean_insert + 3*devi_insert;
	double local_cov;
	int vsize=vbrkpnts.size();
	for(int vi=0;vi<vsize;vi++)
	{
		pos=vbrkpnts[vi].pos;
		//for each candidate site collect features
		//1. left region and right region coverage
		double tmp_lcov, tmp_rcov;
		int region_len=READ_LENGTH/2;
		int region_bases=0;
		for(int i=pos-region_len; i<pos;i++)
		{
			region_bases+=bases[i];
		}
		tmp_lcov=(double)region_bases/(double)region_len;
		
		//cout<<"left coverage: "<<tmp_lcov<<endl;//////////////////////////////////////////////////////////////////////////////////

		int rend=pos+region_len;
		region_bases=0;
		for(int i=pos+1;i<rend;i++)
		{
			region_bases+=bases[i];
		}
		tmp_rcov=(double)region_bases/(double)region_len;

		//cout<<"right coverage: "<<tmp_lcov<<endl;///////////////////////////////////////////////////////////////////////////////////

		if(vbrkpnts[vi].lcov<0.0)
		{//set direction according to coverage, for lclip==rclip situation.
			if(tmp_lcov < tmp_rcov)
				vbrkpnts[vi].bdonor=false;
			else
				vbrkpnts[vi].bdonor=true;
		}

		vbrkpnts[vi].lcov=tmp_lcov;
		vbrkpnts[vi].rcov=tmp_rcov;
		
		//2. fully mapped reads
		//fully_map=base_full_map[pos];
		//vbrkpnts[vi].full_map=base_full_map[pos];
		
		//discordant paired-end reads
		if(vbrkpnts[vi].bdonor==true)
		{//should check right discordant 
			int itemp=0;
			int idx_start=0-discor_range;
			for(int l=0;l<discor_range;l++)
			{
				itemp+=rdiscor[pos-l];
				itemp+=ldiscor[pos-l];
			}
			vbrkpnts[vi].discor=itemp;
		}
		else
		{//should check left discordant  
			int itemp=0;
			for(int l=0;l<discor_range;l++)
			{
				itemp+=ldiscor[pos+l];
				itemp+=rdiscor[pos+l];
			}
			vbrkpnts[vi].discor=itemp;
		}

	}//end of for

	delete[] bases;
	delete[] ldiscor;
	delete[] rdiscor;
	cout<<"step3 finished--------------------------------------------------------------"<<endl;/////////////////////////////////////////////////

	int * clip_pos_order=new int[chrom_len+1]; // give order of each site in chrom_len range, to match vector.
	memset(clip_pos_order,0,chrom_len+1);
	int site_order=0;
	for(int i=0;i<chrom_len;i++)
	{
		if(is_clip_pos[i]==NO_CLIP_POS)
			clip_pos_order[i]=-1;
		else if(is_clip_pos[i]==0)
		{
			clip_pos_order[i]=site_order;
			site_order++;
		}
		else
		{
			if(is_clip_pos[i]>0)
				clip_pos_order[i]=site_order;
			else
				clip_pos_order[i]=site_order-1;
		}
	}

	//output order////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/*ofstream fout_order;
	fout_order.open("clip_order.txt");
	for(int i=0;i<chrom_len;i++)
	{
		if(clip_pos_order[i]!=-1)
		{
			fout_order<<i<<" "<<clip_pos_order[i]<<endl;
		}
	}
	fout_order.close();*/
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	cout<<"step4 finished--------------------------------------------------------------"<<endl;/////////////////////////////////////////////////

	int brkpnt_size=vbrkpnts.size();
	this->nbrkpnt=brkpnt_size;
	this->brkcon=new  BrkConnectNode*[brkpnt_size];
	for(int i=0;i<brkpnt_size;i++)
	{
		brkcon[i]=new BrkConnectNode;
		brkcon[i]->id=i;
		brkcon[i]->pos=vbrkpnts[i].pos;
		brkcon[i]->pnext=NULL;
		brkcon[i]->split_support=0;
	}

	alignSoftClip(fref, fname_clip, is_clip_pos, clip_pos_order);//re-align soft-clipped reads
	delete[] clip_pos_order;
	delete[] is_clip_pos;
	
	cout<<"step5 finished--------------------------------------------------------------"<<endl;/////////////////////////////////////////////////
}


void CandidateSitesCaller::outputFeatures(std::string fout_path)
{
	int vsize=vbrkpnts.size();
	int pos1, pos2;
	int nsplit_map;
	double cov_diff;
	int icanon;
	ofstream fout_junctions;
	fout_junctions.open(fout_path.c_str());

	string fbrkpnt=fout_path+"_brkpnt.txt";
	ofstream fout_brkpnt;
	fout_brkpnt.open(fbrkpnt.c_str());

	vector<JunctionFeatures> vjf;
	vjf.clear();

	string sdirection;//
	for(int i=0;i<vsize;i++)
	{
		double local_cov;
		//output brkpnt information
		if(vbrkpnts[i].bdonor==true) 
		{
			sdirection="-";
			local_cov=vbrkpnts[i].lcov;
		}
		else 
		{
			sdirection="+";
			local_cov=vbrkpnts[i].rcov;
		}
		int nclip=vbrkpnts[i].bdonor==true ? vbrkpnts[i].rclip:vbrkpnts[i].lclip; 
		
		//normalization and output
		fout_brkpnt.setf(ios::fixed);
		fout_brkpnt.precision(4);
		fout_brkpnt<<vbrkpnts[i].pos<<" "<<(double)nclip/local_cov<<" "<<(double)(abs(vbrkpnts[i].lcov-vbrkpnts[i].rcov))/local_cov \
			<<" "<<(double)vbrkpnts[i].softclip_map/local_cov<<" " \
			<<(double)vbrkpnts[i].discor/local_cov<<" "<<sdirection<<endl;
		
		//output junction infomation
		pos1=vbrkpnts[i].pos;
		for(int j=0;j<vbrkpnts[i].vconnct.size();j++)
		{
			int index = vbrkpnts[i].vconnct[j].first;
			pos2 = vbrkpnts[index].pos;
			cov_diff= abs(vbrkpnts[i].lcov-vbrkpnts[i].rcov) + abs(vbrkpnts[index].lcov-vbrkpnts[index].rcov);
			
			if(vbrkpnts[i].vcanon[j]==true)
				icanon=1;
			else
				icanon=0;

			int all_discor=vbrkpnts[i].discor + vbrkpnts[index].discor;
			JunctionFeatures jf;
			if(pos1 > pos2)
			{
				jf.start=pos2; 
				jf.end=pos1;
			}
			else
			{
				jf.start=pos1; 
				jf.end=pos2;
			}

			jf.cov_diff=cov_diff;
			jf.discor=all_discor;
			jf.softclip_map=vbrkpnts[i].vconnct[j].second;
			jf.iscanon=icanon;
			vjf.push_back(jf);
		}//end of for j 
	}//end of for i

	//sort 
	std::sort(vjf.begin(),vjf.end()); 

	//combine and output 
	int vjf_size=vjf.size();
	int last_support, last_start, last_end,  last_iscanon, last_discor;
	double last_covdiff, cur_covdiff;
	int	cur_support, cur_start, cur_end, cur_iscanon, cur_discor;
	last_start=vjf[0].start;
	last_end=vjf[0].end;
	last_support=vjf[0].softclip_map;
	last_iscanon=vjf[0].iscanon;
	last_covdiff=vjf[0].cov_diff;
	last_discor=vjf[0].discor;
	bool bcur=false;
	for(int i=1;i<vjf_size;i++)
	{
		bcur=false;
		cur_start = vjf[i].start;
		cur_end =  vjf[i].end;
		cur_support=vjf[i].softclip_map;
		cur_iscanon=vjf[i].iscanon;
		cur_covdiff=vjf[i].cov_diff;
		cur_discor=vjf[i].discor;
		if((cur_start==last_start) && (cur_end==last_end))
		{
			bcur=true;
			last_support+=cur_support;
		}
		else
		{
			fout_junctions<<last_start<<" "<<last_end<<" "<<last_support<<" "<<last_discor<<" "<<last_covdiff<<" "<<last_iscanon<<endl;
			last_support=cur_support;
		}
		last_start=cur_start;
		last_end=cur_end;
		last_support=cur_support; 
		last_iscanon=cur_iscanon;
		last_discor=cur_discor;
		last_covdiff=cur_covdiff;
	}
	if(bcur==false)
		fout_junctions<<last_start<<" "<<last_end<<" "<<last_support<<" "<<last_discor<<" "<<last_covdiff<<" "<<last_iscanon<<endl;

	fout_junctions.close();
	fout_brkpnt.close();
}

void CandidateSitesCaller::outputGraph()
{
	ofstream fout_graph;
	fout_graph.open(fname_graph.c_str(), ofstream::out);
	//fout_graph<<this->nbrkpnt<<endl;//save total number of breakpoints
	for(int i=0;i<this->nbrkpnt;i++)
	{
		BrkConnectNode* ptmp=brkcon[i];
		while(ptmp!=NULL)
		{
			fout_graph<<ptmp->id<<" ";
			ptmp=ptmp->pnext;
		}
		fout_graph<<endl;
	}
	fout_graph.close();
}

void CandidateSitesCaller::releaseGraph()
{
	//release brkcon 
	for(int i=0;i<this->nbrkpnt;i++)
	{
		BrkConnectNode* ptmp=brkcon[i];
		while(ptmp!=NULL)
		{
			BrkConnectNode* pnode=ptmp->pnext;
			delete ptmp;
			ptmp=pnode;
		}
	}
	delete[] brkcon;
}

void CandidateSitesCaller::loadHardClipRawReads(map<string, string>& lhclip_map, map<string, string>& rhclip_map)
{
	//load left hard-clip
	ifstream fin_lhclip;
	fin_lhclip.open(fname_lhclip_raw.c_str());
	string id, seq;
	while(fin_lhclip>>id>>seq)
	{
		lhclip_map[id]=seq;	
	}
	fin_lhclip.close();

	//load right hard-clip
	ifstream fin_rhclip;
	fin_rhclip.open(fname_rhclip_raw.c_str());
	while(fin_rhclip>>id>>seq)
	{
		rhclip_map[id]=seq;	
	}
	fin_rhclip.close();

}

void CandidateSitesCaller::alignSoftClip(std::string fref, std::string freads, int*& is_clip_pos, int*& clip_pos_order)
{
	ifstream fin_reads;
	fin_reads.open(freads.c_str());

	ReferenceSeq rseq(fref);
	string chrom_name,scigar, seq, qname;
	int clip_type, map_pos, pnext, flag;
	int len1, len2, clip_pos1, clip_pos2;
	Alignment almt;
	
	string sfront="";
	string stail="";
	int cnt_reads=0;
	
	map<string, string> lhclip_map;
	lhclip_map.clear();
	map<string, string> rhclip_map;
	rhclip_map.clear();
	loadHardClipRawReads(lhclip_map, rhclip_map);

	//ofstream ftest_hclip;///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//ftest_hclip.open("test_hard_clip.txt");/////////////////////////////////////////////////////////////////////////////////////////////
	ofstream funmapped_sclip;
	funmapped_sclip.open("unmapped_soft_clip_reads.txt");

	while(fin_reads>>clip_type>>qname>>flag>>chrom_name>>map_pos>>scigar>>seq>>pnext \
						>>len1>>len2>>clip_pos1>>clip_pos2)
	{
		//output hint information
		if((cnt_reads%10000)==0) 
		{
			cout<<cnt_reads<<" reads has been re-aligned!!!"<<endl;
		}
		cnt_reads++;

		//replace to original or bypass hard-clip 
		if(clip_type==READ_TYPE_LEFT_HARDCLIP || clip_type==READ_TYPE_RIGHT_HARDCLIP || clip_type==READ_TYPE_BOTH_HARDCLIP)
		{
			//first check whether can replace
			bool breplce=false;
			if((flag&0x40) != 0)
			{//is first in pair
				map<string, string>::iterator iter = lhclip_map.find(qname);
				if(iter != lhclip_map.end())
				{//find one
					seq = iter->second;
					breplce=true;
					//ftest_hclip<<qname<<" "<<seq<<endl;/////////////////////////////////////////////////////////////////////////////////////////////
				}
			}

			if((flag&0x80) != 0)
			{//is second in pair
				map<string, string>::iterator iter = rhclip_map.find(qname);
				if(iter != rhclip_map.end())
				{//find one
					seq = iter->second;
					breplce=true;
					//ftest_hclip<<qname<<" "<<seq<<endl;/////////////////////////////////////////////////////////////////////////////////////////////
				}
			}
			//if not, then bypass
			if(breplce==false) continue;
		}

		vector<pair<string,int> > vcigar;
		vcigar.clear();
		almt.str2Cigar(scigar,vcigar);
		
		int pos_front, len_front;//position on read and length of clip
		int pos_tail, len_tail;
		int seg_len=seq.length();
		//cout<<"segment info:"<< pos_front<<" "<<len_front<<" "<<pos_tail<<" "<<len_tail<<endl;/////////////////////////////////////////////////////////////
		getCigarPosLen(vcigar,seg_len,pos_front,len_front,pos_tail,len_tail);//only consider front and tail maybe not that precisely/////////////////
		vector<pair<string,int> >().swap(vcigar);//release vcigar; 

		int vi;
		bool bhit=false;
		if(pos_front!=-1 && is_clip_pos[clip_pos1] != NO_CLIP_POS)
		{// front segment clipped 
			if(len_front<CUTOFF_LA_SGMT) continue;
			sfront=seq.substr(pos_front,len_front);
			vi=clip_pos_order[clip_pos1];
			bhit=searchMapSite(sfront,rseq, chrom_name, vi, pnext, true);
			if(bhit==false)
			{
				funmapped_sclip<<qname<<" "<<flag<<" "<<scigar<<" "<<seq<<" "<<clip_pos1<<" "<<clip_pos2<<endl;
				if(sfront=="GCAATTTCTCAA") cout<<qname<<" "<<flag<<" "<<scigar<<" "<<seq<<" "<<clip_pos1<<" "<<clip_pos2<<endl;
			}
		}

		if(pos_tail!=-1 && is_clip_pos[clip_pos2] != NO_CLIP_POS)
		{//tail segment clipped 
			if(len_tail<CUTOFF_LA_SGMT) continue;
			stail=seq.substr(pos_tail,len_tail);
			vi=clip_pos_order[clip_pos2];
			bhit=searchMapSite(stail,rseq, chrom_name, vi, pnext, false);
			if(bhit==false)
			{
				funmapped_sclip<<qname<<" "<<flag<<" "<<scigar<<" "<<seq<<" "<<clip_pos1<<" "<<clip_pos2<<endl;
				if(stail=="GCAATTTCTCAA") cout<<qname<<" "<<flag<<" "<<scigar<<" "<<seq<<" "<<clip_pos1<<" "<<clip_pos2<<endl;
			}
		}
	}

	funmapped_sclip.close();
	map<string, string>().swap(lhclip_map);
	map<string, string>().swap(rhclip_map);

	//ftest_hclip.close();////////////////////////////////////////////////////////////////////////////////////////////////////////////
	fin_reads.close();
}

bool CandidateSitesCaller::searchMapSite(std::string segmnt, ReferenceSeq& refseq, std::string chrom_name, int vi, int pnext, bool bfront)
{
	//cout<<"search segment: "<<segmnt<<" "<<vi<<endl;////////////////////////////////////////////////////////////////////////////////////////////// 
	int clip_pos=vbrkpnts[vi].pos;
	bool isdonor= vbrkpnts[vi].bdonor;
	//cout<<"clip pos: "<<clip_pos<<endl;/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	int lbound=0;
	int rbound=0;
	std::string sub_ref;
	int ref_pos, ref_length;

	if((isdonor==true && bfront==true) || (isdonor==false && bfront==false))
	{
		//this read is wrongly clipped here, 
		//some process here ...
		cout<<"wrong read: "<<endl;/////////////////////////////////////////////////////////////////////////////////////////////////
		return false;
	}

	bool bhit=false;
	if(bfront==true)
	{//front segment clipped, and brkpnt is acceptor site
		//cout<<"Front segment clipped"<<endl;////////////////////////////////////////////////////////////////////////////////////////////////
		if(pnext<clip_pos)
		{//use pnext as bound
			lbound=pnext;
		}
		else
		{
			lbound=clip_pos-MAX_INTRON_SIZE;
		}

		int index=vi-1;
		//cout<<"Boundaries: "<<lbound<<" "<<index<<endl;/////////////////////////////////////////////////////////////////////////////////
		while(index >=0 && vbrkpnts[index].pos > lbound)
		{
			//check brkpnt direction
			if(vbrkpnts[index].bdonor==false)
			{// current is acceptor site. Bypass to find an donor site.
				index--; 
				continue;
			}
			
			if(vbrkpnts[index].pos==0) cout<<"index "<<index<<" "<<vbrkpnts[index].lclip<<" "<<vbrkpnts[index].rclip<<endl;

			//get reference substring 
			//cout<<"Get front ref-------------------------------------------------"<<endl;///////////////////////////////////////////////
			ref_pos=vbrkpnts[index].pos-BASE_SLACK;
			ref_length=segmnt.length() + 2*BASE_SLACK;
			//cout<<"Front, try to align with "<<ref_pos<<" "<<ref_length;//////////////////////////////////////////////////////////////////
			refseq.getSubRef(chrom_name,ref_pos,ref_length, sub_ref);
			//cout<<" "<<sub_ref<<endl;/////////////////////////////////////////////////////////////////////////////////////////////////////

			if(segmnt=="GCAATTTCTCAA") cout<<"front_clip "<<lbound<<" "<<clip_pos<<" "<<vbrkpnts[index].pos<<" "<<sub_ref<<endl;///////////////////////////////////////////////////////////////////////////////////////
			//check whether can be fully mapped 
			bhit=alignSegmt(sub_ref, segmnt);
			if(segmnt=="GCAATTTCTCAA" && bhit==true) cout<<segmnt<<" hit!!!!!!"<<endl;
			//cout<<"step3.0=-------------------------------------------------------------"<<endl;///////////////////////////////////////////////
			if(bhit==true)
			{//find a hit 
				//some process
				//need to check the mapping position ???????????????????????????

				//////information of the graph///////////////////////////////////////////////////////////////////////////////////
				//save to this candidate site a record 
				bool bexist=false;
				BrkConnectNode* pnode=this->brkcon[vi];
				while(pnode!=NULL)
				{
					if(pnode->id==index)
					{
						pnode->split_support++;
						bexist=true;
						break;
					}
					pnode=pnode->pnext;
				}
				
				if(bexist==false)
				{
					BrkConnectNode* pnew= new BrkConnectNode;
					pnew->id=index;
					pnew->pos=vbrkpnts[index].pos;
					pnew->split_support=1;
					pnew->pnext=this->brkcon[vi]->pnext;
					this->brkcon[vi]->pnext=pnew;
				}
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////

				bexist=false;
				for(int k=0; k<vbrkpnts[vi].vconnct.size(); k++) /////////need better algorithm to search, like binary search
				{
					if(vbrkpnts[vi].vconnct[k].first==index)
					{
						vbrkpnts[vi].vconnct[k].second++;
						bexist=true;
						break;
					}
				}
				
				if(bexist==false)
				{				
					vbrkpnts[vi].vconnct.push_back(std::make_pair(index,1));				
					//check whether is canonical from reference
					string slref, srref;
					refseq.getSubRef(chrom_name,clip_pos-CANON_SLACK,2+2*CANON_SLACK, srref);
					size_t offset;
					bool is_acceptor=this->isAcceptorSite(srref,offset);
					refseq.getSubRef(chrom_name,vbrkpnts[vi].pos-CANON_SLACK,2+2*CANON_SLACK, slref);
					bool is_donor=this->isDonorSite(slref,offset);
					if(is_acceptor==true && is_donor==true) 
						vbrkpnts[vi].vcanon.push_back(true);
					else 
						vbrkpnts[vi].vcanon.push_back(false);
				}

				vbrkpnts[index].softclip_map++;
				break;//find a hit
			}
			else
			{
				index--;
			}
		}//end of while
	}
	else
	{//tail segment clipped //donor site 
		//cout<<"Tail segment clipped"<<endl;//////////////////////////////////////////////////////////////////////////////////////
		if(pnext > clip_pos)
		{
			rbound=pnext;
		}
		else
		{
			rbound=clip_pos + MAX_INTRON_SIZE;
		}
		int index=vi+1;
		//cout<<"Boundaries: "<<rbound<<" "<<index<<endl;/////////////////////////////////////////////////////////////////////////////////
		while((index< vbrkpnts.size()) && (vbrkpnts[index].pos < rbound))
		{
			if(vbrkpnts[index].bdonor==true) //only check donor site, if not, bypass
			{//acceptor site
				index++; 
				continue;
			}
			
			ref_pos=vbrkpnts[index].pos-BASE_SLACK;
			ref_length=segmnt.length() + 2*BASE_SLACK;
			//cout<<"Tail, try to align with "<<ref_pos<<" "<<ref_length;///////////////////////////////////////////////
			refseq.getSubRef(chrom_name,ref_pos,ref_length, sub_ref);

			//cout<<" "<<sub_ref<<endl;////////////////////////////////////////////////////////////////////////////////////
			if(segmnt=="GCAATTTCTCAA") cout<<"tail_clip "<<rbound<<" "<<clip_pos<<" "<<vbrkpnts[index].pos<<" "<<sub_ref<<endl;///////////////////////////////////////////////////////////////////////////////////////

			//check whether can be fully mapped check
			bhit=alignSegmt(sub_ref, segmnt);

			if(segmnt=="GCAATTTCTCAA" && bhit==true) cout<<segmnt<<" hit!!!!!!"<<endl;

			if(bhit==true)
			{//find a hit
				//some process 
				//need to check the mapping position ???????????????????????????

				//////information of the graph///////////////////////////////////////////////////////////////////////////////////
				//save to this candidate site a record
				bool bexist=false;
				BrkConnectNode* pnode=this->brkcon[vi];
				while(pnode!=NULL)
				{
					if(pnode->id==index)
					{
						pnode->split_support++;
						bexist=true;
						break;
					}
					pnode=pnode->pnext;
				}
				
				if(bexist==false)
				{
					BrkConnectNode* pnew= new BrkConnectNode;
					pnew->id=index;
					pnew->pos=vbrkpnts[index].pos;
					pnew->split_support=1;
					pnew->pnext=this->brkcon[vi]->pnext;
					this->brkcon[vi]->pnext=pnew;
				}
				///////////////////////////////////////////////////////////////////////////////////////////////////////////

				bexist=false;
				for(int k=0;k<vbrkpnts[vi].vconnct.size(); k++)
				{
					if(vbrkpnts[vi].vconnct[k].first==index)
					{
						vbrkpnts[vi].vconnct[k].second++;
						bexist=true;
						break;
					}
				}
				
				if(bexist==false) 
				{
					vbrkpnts[vi].vconnct.push_back(std::make_pair(index,1));
					
					//check whether is canonical from reference
					string slref, srref;				
					refseq.getSubRef(chrom_name,clip_pos-CANON_SLACK,2+2*CANON_SLACK, slref);
					refseq.getSubRef(chrom_name,vbrkpnts[vi].pos-CANON_SLACK,2+2*CANON_SLACK, srref);
					size_t offset;
					bool is_acceptor=this->isAcceptorSite(srref,offset);
					bool is_donor=this->isDonorSite(slref,offset);
					if(is_acceptor==true && is_donor==true) vbrkpnts[vi].vcanon.push_back(true);
					else vbrkpnts[vi].vcanon.push_back(false);
				}
				vbrkpnts[index].softclip_map++;
				break;
			}
			else
			{
				index++;
			}
		}
	}//end of else 
	//cout<<"End of Search"<<endl;/////////////////////////////////////////////////////////////////////////////////////////// 
	return bhit;
}


bool CandidateSitesCaller::alignSegmt(string sub_ref, string str_clip)
{
	LocalAlignment la;
	la.setRef(sub_ref);
	la.setSgmt(str_clip);
	int opt_ref_start, opt_ref_end, opt_sgmt_start, opt_sgmt_end;
	la.optAlign(opt_ref_start, opt_ref_end, opt_sgmt_start,opt_sgmt_end);

	int len_mapped=opt_sgmt_end-opt_sgmt_start+1;
	double hit=(double)len_mapped/(double)str_clip.length();
	//if(str_clip=="GCAATTTCTCAA") cout<<hit<<" "<<len_mapped<<endl;///////////////////////////////////////////////////////////////////////////////
	if(hit >= HIT_LOCAL_ALIGNMENT)
		return true;
	else
		return false;
}


/*
Description:
	For soft-clip or hard-clip reads, get position and length of clip-part.
Input:
	cigar: cigar of reads.
Output:
	pos_front: when clip happens at front part, o/w -1.
	len_front: length of the front clip part, o/w -1.
	pos_tail: when clip happens at tail part, o/w -1.
	len_tail: length of the tail clip part, o/w -1.
*/
void CandidateSitesCaller::getCigarPosLen(std::vector< std::pair<std::string,int> >& cigar, int seg_len, int& pos_front, int& len_front, int& pos_tail, int& len_tail)
{	
	int size=cigar.size();

	int cnt=0;
	pos_front=-1;
	len_front=-1;
	pos_tail=-1;
	len_tail=-1;
	
	if(cigar[0].first=="S" || cigar[0].first=="H")
	{
		pos_front=0;
		len_front=cigar[0].second;
	}

	if(size>1 && (cigar[size-1].first=="S" || cigar[size-1].first=="H"))
	{		
		len_tail=cigar[size-1].second;
		pos_tail=seg_len-len_tail;
	}
}


void CandidateSitesCaller::setRefPath(std::string fref)
{
	this->fref=fref;
}

//another pair is GC-AG
//left of a intron is the donor splice site
bool CandidateSitesCaller::isDonorSite(string str, std::size_t& offset)
{
	offset=str.find("GT");
	if(offset==std::string::npos)
		return false;
	else 
		return true;
}

//right of intron is the acceptor splice site 
bool CandidateSitesCaller::isAcceptorSite(string str, std::size_t& offset)
{
	offset=str.find("AG");
	if(offset==std::string::npos)
		return false;
	else 
		return true;
}

void CandidateSitesCaller::setInsertSize(double mean, double devi)
{
	this->mean_insert=mean;
	this->devi_insert=devi;
}
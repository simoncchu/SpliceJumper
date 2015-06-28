#include"Coverage.h"
#include"bam_parse.h"
#include"Alignment.h"
#include"public_parameters.h"
#include<iostream>
#include<cstring>
#include<fstream>
#include"CandidateSitesCaller.h"

Coverage::Coverage()
{
	bam_path="";
}

Coverage::Coverage(std::string bam_name, std::string fai_path)
{
	this->bam_path=bam_name;
}

void Coverage::initFiles()
{
	ofstream fout_fm;
	fout_fm.open(fname_fm.c_str(),ofstream::out);
	ofstream fout_clip;
	fout_clip.open(fname_clip.c_str(),ofstream::out);
	ofstream fout_um;
	fout_um.open(fname_um.c_str(),ofstream::out);
	ofstream fout_others;
	fout_others.open(fname_others.c_str(),ofstream::out);
	ofstream fout_cov;
	fout_cov.open(fname_cov.c_str(),ofstream::out);
	ofstream fout_brkpnts;
	fout_brkpnts.open(fname_brkpnts.c_str(), ofstream::out);
	ofstream fout_lhclip_id;
	fout_lhclip_id.open(fname_left_hclip.c_str(),ofstream::out);
	ofstream fout_rhclip_id;
	fout_rhclip_id.open(fname_right_hclip.c_str(),ofstream::out);

	fout_fm.close();
	fout_clip.close();
	fout_um.close();
	fout_others.close();
	fout_cov.close();
	fout_brkpnts.close();
	fout_lhclip_id.close();
	fout_rhclip_id.close();
}

//--------public functions-------------------------------------------------------
/*
Description:
	calc the coverage.
Note: make sure the bam file is sorted by chromosome, mapping position.
*/
void Coverage::pileup(int cid, BamParse& bp)
{
	string chrom_name;
	int chrom_len=-1;

	ofstream fout_fm;
	fout_fm.open(fname_fm.c_str(),ofstream::out|ofstream::app);
	ofstream fout_clip;
	fout_clip.open(fname_clip.c_str(),ofstream::out|ofstream::app);
	ofstream fout_um;
	fout_um.open(fname_um.c_str(),ofstream::out|ofstream::app);
	ofstream fout_others;
	fout_others.open(fname_others.c_str(),ofstream::out|ofstream::app);
	ofstream fout_cov;
	fout_cov.open(fname_cov.c_str(),ofstream::out|ofstream::app);
	ofstream fout_brkpnts;
	fout_brkpnts.open(fname_brkpnts.c_str(), ofstream::out|ofstream::app);
	
	ofstream fout_lhclip_id;
	fout_lhclip_id.open(fname_left_hclip.c_str(),ofstream::out|ofstream::app);
	ofstream fout_rhclip_id;
	fout_rhclip_id.open(fname_right_hclip.c_str(),ofstream::out|ofstream::app);

	int window=20000000; 
	int start=0;
	int end=window;
	
	bp.getChromNameLength(cid, chrom_name, chrom_len);
	if(chrom_len<0)
	{
		cout<<"Can not load chromosome: "<<chrom_name<<endl;
		return;
	}

	int* bases=new int[chrom_len+1];//save the coverage of each chromosome
	//memset(bases,0,chrom_len+1);
	int* is_clip_pos=new int[chrom_len+1]; //record whether this position is an candidate position, combined ones record the distance to the record one.
	//memset(is_clip_pos,NO_CLIP_POS,chrom_len+1);
	int* ldiscor=new int[chrom_len+1];
	//memset(ldiscor,0,chrom_len+1);
	int* rdiscor=new int[chrom_len+1];
	//memset(rdiscor,0,chrom_len+1);
	Brkpnt* brkpnts=new Brkpnt[chrom_len+1];
	for(int i=0;i<chrom_len;i++)
	{
		bases[i]=0;
		is_clip_pos[i]=NO_CLIP_POS;
		ldiscor[i]=0;
		rdiscor[i]=0;
	}

	cout<<chrom_name<<" "<<chrom_len<<endl;//////////////////////////////////////////////////////////
	bool bswitch=false;
	while(end<=chrom_len)
	{
		cout<<"start: "<<start<<" end: "<<end<<endl;/////////////////////////////////////////////////
		pileupRegion(cid, start,end, bp, bases, brkpnts, is_clip_pos, ldiscor, rdiscor,\
			fout_fm, fout_clip, fout_um,fout_others,fout_cov, fout_lhclip_id, fout_rhclip_id, chrom_name);
		start=end+1;
		end+=window;
		if(end>chrom_len)
		{
			if(bswitch==true) break;
			else
			{
				bswitch=true;
				end=chrom_len;
			}
		}
	}

	cout<<"Check discordant 4483852: "<<rdiscor[4483852]<<endl;

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//ofstream fout_test_discor;///////////////////////////////////////////////////////////////////////////////////////////////////
	//fout_test_discor.open("discor_test.txt");
	//for(int i=0;i<chrom_len;i++)
	//{
	//	if(rdiscor[i]!=0)
	//		fout_test_discor<<i<<" "<<rdiscor[i]<<endl;
	//}
	//fout_test_discor.close();///////////////////////////////////////////////////////////////////////////////////////////////////

	//output the coverage information
	for(int i=0;i<chrom_len;i++)
	{
		if(bases[i]>0)
			fout_cov<<chrom_name<<"\t"<<i<<"\t"<<bases[i]<<endl;
	}
		
	vector< Brkpnt > vcandidates;
	vcandidates.clear();
	
	callCandidateByCov(brkpnts, bases, is_clip_pos, ldiscor, rdiscor, chrom_len);
	
	//output the brkpnts information 
	for(int i=0;i<chrom_len;i++)
	{
		int sum_clip=brkpnts[i].lclip+brkpnts[i].rclip;
		if(sum_clip>=MIN_CLIP_READS) //at least 1 reads clipped at the same position.
		{//clipped sites
			vcandidates.push_back(brkpnts[i]);
		}
		else if(sum_clip < 0)
		{//call out by coverage, no clip happen there.
			vcandidates.push_back(brkpnts[i]);
		}
		else
		{
			is_clip_pos[i]=NO_CLIP_POS;
		}
	}
	delete[] brkpnts;

	//combine brkpnts close to each other.
	int vsize=vcandidates.size();
	int last_pos=vcandidates[0].pos;
	int last_lclip=vcandidates[0].lclip;
	int last_rclip=vcandidates[0].rclip;
	bool last_bdonor, cur_bdonor;
	if(last_lclip>last_rclip) last_bdonor=false;
	else last_bdonor=true;

	int cur_pos, cur_lclip, cur_rclip;
	//int pos=last_pos, split=last_split;
	for(int i=1;i<vsize;i++)
	{
		cur_pos=vcandidates[i].pos;
		cur_lclip=vcandidates[i].lclip;
		cur_rclip=vcandidates[i].rclip;
		if(cur_lclip>cur_rclip) cur_bdonor=false;
		else cur_bdonor=true;

		if((cur_pos-last_pos > BASE_SLACK) || (cur_bdonor!=last_bdonor))
		{//no need to combine
			string sdir;
			if(last_lclip>last_rclip)
			{
				sdir="+";
				if(last_rclip==-1) last_rclip=0; //this happen when call callCandidateByCov() function
			}
			else if(last_lclip<last_rclip)
			{
				sdir="-";
				if(last_lclip==-1) last_lclip=0; //this happen when call callCandidateByCov() function 
			}
			else
			{
				sdir="*";
			}

			fout_brkpnts<<chrom_name<<" "<<last_pos<<" "<<last_lclip<<" "<< last_rclip<<" "<<sdir<<endl;
			last_pos=cur_pos;
			last_lclip=cur_lclip;
			last_rclip=cur_rclip;
			last_bdonor=cur_bdonor;
		}
		else
		{//need to combine
			int max_pos=last_pos;
			int max_clip=0;
			int max_i=i-1;

			int tmp_i=i; 

			if(last_bdonor==true)
			{
				max_clip=last_rclip;
			}
			else
			{
				max_clip=last_lclip;	
			}

			while(tmp_i<vsize)
			{//in the range and same direction
				cur_pos=vcandidates[tmp_i].pos;
				cur_lclip=vcandidates[tmp_i].lclip;
				cur_rclip=vcandidates[tmp_i].rclip;
				if(cur_lclip>cur_rclip) cur_bdonor=false;
				else cur_bdonor=true;

				if((cur_pos-last_pos > BASE_SLACK) || (cur_bdonor!=last_bdonor)) 
				{//
					for(int k=i-1;k<tmp_i;k++)
					{
						is_clip_pos[vcandidates[k].pos]=max_pos-vcandidates[k].pos;
					}
					///////////////////////////////////////////////////////////////////////////////////////////////
					string sdir;
					if(vcandidates[max_i].lclip > vcandidates[max_i].rclip)
					{
						sdir="+";
						if(vcandidates[max_i].rclip==-1) vcandidates[max_i].rclip=0; //this happen when call callCandidateByCov() function
					}
					else if(vcandidates[max_i].lclip < vcandidates[max_i].rclip)
					{
						sdir="-";
						if(vcandidates[max_i].lclip==-1) vcandidates[max_i].lclip=0; //this happen when call callCandidateByCov() function 
					}
					else
					{// no direction
						sdir="*";
					}
					//////////////////////////////////////////////////////////////////////////////////////////////////////////
					fout_brkpnts<<chrom_name<<" "<<vcandidates[max_i].pos<<" "<<vcandidates[max_i].lclip<<" "<<vcandidates[max_i].rclip<<" "<<sdir<<endl;
					
					last_pos=cur_pos;
					last_lclip=cur_lclip;
					last_rclip=cur_rclip;
					last_bdonor=cur_bdonor;
					break;
				}
				else
				{//combine 
					if(cur_bdonor==true)
					{
						if(cur_rclip > max_clip)
						{
							max_clip=cur_rclip;
							max_pos=vcandidates[tmp_i].pos;
							max_i=tmp_i;
						}
					}
					else
					{
						if(cur_lclip > max_clip)
						{
							max_clip=cur_lclip;
							max_pos=vcandidates[tmp_i].pos;
							max_i=tmp_i;
						}
					}
				}
				last_pos=cur_pos;
				last_lclip=cur_lclip;
				last_rclip=cur_rclip;
				last_bdonor=cur_bdonor;
				tmp_i++;
			}//end of while 
			i=tmp_i;
		}
	}
	
	string sdir;
	if(last_lclip>last_rclip) sdir="+";
	else if(last_lclip < last_rclip)sdir="-";
	else sdir="*";
	fout_brkpnts<<chrom_name<<" "<<last_pos<<" "<<last_lclip<<" "<< last_rclip<<" "<<sdir<<endl;

	saveBrkpntInfo(chrom_len,is_clip_pos,ldiscor,rdiscor);//save brkpnt information into file
	
	delete[] bases;
	delete[] is_clip_pos;
	delete[] ldiscor;
	delete[] rdiscor;

	fout_fm.close();
	fout_clip.close();
	fout_um.close();
	fout_others.close();
	fout_cov.close();
	fout_brkpnts.close();
	//fout_hclip.close();
	fout_lhclip_id.close();
	fout_rhclip_id.close();
	
}


void Coverage::saveBrkpntInfo(int chrom_len, int*& is_clip_pos, int*& ldiscor, int*& rdiscor)
{
	ofstream fout_icp;
	fout_icp.open(fname_icp.c_str());
	ofstream fout_ld;
	fout_ld.open(fname_ldiscor.c_str());
	ofstream fout_rd;
	fout_rd.open(fname_rdiscor.c_str());
	for(int i=0;i<chrom_len;i++)
	{
		if(is_clip_pos[i]!=NO_CLIP_POS)
			fout_icp<<i<<" "<<is_clip_pos[i]<<endl;
		if(ldiscor[i]!=0)
			fout_ld<<i<<" "<<ldiscor[i]<<endl;
		if(rdiscor[i]!=0)
			fout_rd<<i<<" "<<rdiscor[i]<<endl;
	}
	fout_rd.close();
	fout_ld.close();
	fout_icp.close();
}


void Coverage::pileupRegion(int cid, int start, int end, BamParse& bp, int*& bases, Brkpnt*& brkpnts, int*& is_clip_pos, int*& ldiscor, int*& rdiscor, \
	ofstream& fout_fm, ofstream& fout_clip, ofstream& fout_um,ofstream& fout_others, ofstream& fout_cov, \
	ofstream& fout_lclip, ofstream& fout_rclip, string& chrom_name)
{
	bp.clearAll();//clear all the records 
	bool bsignal=bp.parseAlignment(cid,start,cid,end);//read the whole reads if this chromosome into memory
	if(bsignal==false)
	{
		std::cout<<"Cannot parse bam file "<<this->bam_path<<std::endl;
		return;
	}

	int discor_range=mean_insert + 3*devi_insert;
	int max_range=discor_range + MAX_INTRON_SIZE;
	int vbamsize=bp.bam_aln_records.size();//number of reads
	for(int i=0;i<vbamsize;i++)
	{		
		int map_pos=bp.bam_aln_records[i]->pos;
		if(map_pos<start || map_pos>end) 
		{
			//cout<<"Out of range reads: "<<map_pos<<endl;//////////////////////////////
			continue;
		}

		Alignment alnmt;
		alnmt.setBar(bp.bam_aln_records[i]);
			
		std::string cigar=alnmt.getCigar();
		string seq=bp.bam_aln_records[i]->seq;
		string qname=bp.bam_aln_records[i]->qName;
		int rid=bp.bam_aln_records[i]->rID;
		int rnext_id=bp.bam_aln_records[i]->rNextID;// ref id of mate-read 
		int pnext=bp.bam_aln_records[i]->pNext;//mate read mapping position.
		int flag=bp.bam_aln_records[i]->flag;
		int read_type;
		//first check whether is a qualified read. 
		if((alnmt.isDuplicate()==true) || (alnmt.isPrimaryAlign() == false) || (alnmt.passQualityCK()==false))
		{
			read_type=READ_TYPE_OTHER;
			//write into others 
			fout_others<<read_type<<" "<<qname<<" "<<flag<<" "<<chrom_name<<" "<<map_pos<<" "<<cigar<<" "<<seq<<" "<<pnext<<endl; 
			continue;
		}
		//then check reads type 
		int aln_type=alnmt.getReadType();
		if(aln_type==READ_TYPE_FULLMAP)
		{
			for(int i=0;i<READ_LENGTH;i++)
			{
				bases[map_pos+i]++;
				//full_map[map_pos+i]++;
			}
			//write into file 
			fout_fm<<aln_type<<" "<<qname<<" "<<flag<<" "<<chrom_name<<" "<<map_pos<<" "<<cigar<<" "<<seq<<" "<<pnext<<endl; 
		}
		else if(aln_type==READ_TYPE_CLIP)
		{
			int pos1,len1,pos2,len2, clip_pos1, clip_pos2;
			int clip_type=alnmt.getClipType(pos1,len1,pos2,len2);
			if(clip_type==READ_TYPE_LEFT_SOFTCLIP || clip_type==READ_TYPE_LEFT_HARDCLIP)
			{
				int j=0;
				for(int i=len1;i<READ_LENGTH;i++)
				{
					bases[map_pos+j]++;
					j++;
				}
				clip_pos1=map_pos;
				clip_pos2=-1;
			}
			else if(clip_type==READ_TYPE_RIGHT_SOFTCLIP || clip_type==READ_TYPE_RIGHT_HARDCLIP)
			{
				int len_mapped = READ_LENGTH - len2;
				for(int i=0;i<len_mapped;i++)
				{
					bases[map_pos+i]++;
				}
				clip_pos2 = map_pos + len_mapped;
				clip_pos1=-1;
			}
			else if(clip_type==READ_TYPE_BOTH_SOFTCLIP || clip_type==READ_TYPE_BOTH_HARDCLIP)
			{
				int j=0;
				for(int i=len1;i<READ_LENGTH-len2;i++)
				{
					bases[map_pos+j]++;
					j++;
				}
				clip_pos1=map_pos;
				clip_pos2=map_pos + (READ_LENGTH - len1 - len2); 
			}
						
			if(clip_pos1!=-1)
			{//left soft-clip
				brkpnts[clip_pos1].pos=clip_pos1;
				brkpnts[clip_pos1].lclip++;	
				is_clip_pos[clip_pos1]=0;
			}
			if(clip_pos2!=-1)
			{// right soft-clip
				brkpnts[clip_pos2].pos=clip_pos2;
				brkpnts[clip_pos2].rclip++;
				is_clip_pos[clip_pos2]=0;
			}
			if(clip_type==READ_TYPE_LEFT_HARDCLIP || clip_type==READ_TYPE_RIGHT_HARDCLIP || clip_type==READ_TYPE_BOTH_HARDCLIP)
			{
				//fout_hclip<<clip_type<<" "<<qname<<" "<<flag<<" "<<chrom_name<<" "<<map_pos<<" "<<cigar<<" "<<seq<<" "<<pnext \
						//<<" "<<len1<<" "<<len2<<" "<<clip_pos1<<" "<<clip_pos2<<endl;
				if(alnmt.isFirstInPair()==true)
				{//output reads ids of "first in pair" reads
					fout_lclip<<qname<<endl;
				}

				if(alnmt.isSecondInPair()==true)
				{//output reads ids of "second in pair" reads
					fout_rclip<<qname<<endl;
				}

			}
			//else 
			fout_clip<<clip_type<<" "<<qname<<" "<<flag<<" "<<chrom_name<<" "<<map_pos<<" "<<cigar<<" "<<seq<<" "<<pnext \
					<<" "<<len1<<" "<<len2<<" "<<clip_pos1<<" "<<clip_pos2<<endl; 
			
		}
		else if(aln_type==READ_TYPE_UNMAP)
		{//unmapped reads
			//check mate is mapped or not 
			int type;
			if(alnmt.isMateMapped()==true)
			{
				type=READ_PAIR_MAP_TYPE_01;
			}
			else
			{
				type=READ_PAIR_MAP_TYPE_00;
			}
			fout_um<<type<<" "<<qname<<" "<<flag<<" "<<chrom_name<<" "<<map_pos<<" "<<cigar<<" "<<seq<<" "<<pnext<<endl;//
			
		}
		else
		{//other 
			fout_others<<qname<<" "<<flag<<" "<<chrom_name<<" "<<map_pos<<" "<<cigar<<" "<<seq<<" "<<pnext<<endl;// other reads 
		}

		//check whether discordant or not
		int cor_status=isDiscordant(rid, map_pos,rnext_id, pnext, discor_range, max_range);
		
		if(cor_status==LDISCORDANT)
		{
			ldiscor[map_pos]++;
		}
		else if(cor_status==RDISCORDANT)
		{
			rdiscor[map_pos]++;
		}
		
	}//end of for

	bp.clearAll();
}

/*
Description:
Call candidate sites by coverage change and discordant paired-end reads.
*/
void Coverage::callCandidateByCov(Brkpnt*& brkpnts, int*& cov, int*& is_clip_pos, int*& ldiscor, int*& rdiscor, int chrom_len)
{
	int last_cov=cov[0];
	int cur_cov=0;
	for(int i=1;i<chrom_len;i++)
	{
		cur_cov=cov[i];

		//check whether is a clip site, if yes, then bypass
		if(brkpnts[i].lclip>0 || brkpnts[i].rclip>0)
		{
			last_cov=cur_cov;
			continue;
		}

		//int range=this->mean_insert + 3*this->devi_insert;
		int range=READ_LENGTH/2;
		int cnt=0;
		bool bdonor=true;
		if(last_cov==0 && cur_cov>0) 
		{
			//check whether have discordant paired-end reads support
			for(int j=0;j<range;j++)
			{
				cnt+=ldiscor[i+j];
				cnt+=rdiscor[i+j];
			}
		}
		else if(last_cov>0 && cur_cov==0)
		{
			for(int j=0;j<range;j++)
			{
				cnt+=ldiscor[i-j];
				cnt+=rdiscor[i-j];
			}
			bdonor=false;
		}
		
		if(cnt>0)
		{
			//check whether already exist
			if(brkpnts[i].pos==0)
			{
				brkpnts[i].pos=i;
				is_clip_pos[i]=0;
				if(bdonor==true)
					brkpnts[i].rclip=-1;
				else
					brkpnts[i].lclip=-1;
			}
		}
		
		last_cov=cur_cov;
	}
}

int Coverage::isDiscordant(int rid, int map_pos, int rnext_id, int pnext, int discor_range, int max_range)
{
	if(rid != rnext_id) return WRONG_PAIR;

	int dist=abs(map_pos-pnext);
	if((dist > discor_range) && (dist< max_range))
	{//discordant happen
		if(map_pos> pnext) return LDISCORDANT;
		else return RDISCORDANT;
	}
	else if(dist<discor_range)
	{//concordant happen
		return CONCORDANT;
	}
	else
	{//wrong mapping
		return WRONG_PAIR;		
	}
}


//-------------------public member functions-----------------------------------
std::string Coverage::getBamPath()
{
	return this->bam_path;
}

void Coverage::setBamPath(std::string bam_name)
{
	this->bam_path=bam_name;
}

//std::string Coverage::getFaiPath()
//{
//	return this->fai_path;
//}
//
//void Coverage::setFaiPath(std::string fai_path)
//{
//	this->fai_path=fai_path;
//}


void Coverage::setInsertSize(double mean, double devi)
{
	this->mean_insert=mean;
	this->devi_insert=devi;
}
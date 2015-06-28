#include"ReferenceSeq.h"

ReferenceSeq::ReferenceSeq()
{
	
}

ReferenceSeq::ReferenceSeq(std::string path)
{
	this->path=path;
	fp.setPath(path);
}

void ReferenceSeq::getSubRef(std::string chrom_name,int pos, int length, std::string& sub_ref)//return sub-reference 
{
	sub_ref = this->fp.parseFasta(chrom_name,pos,length);	
}

void ReferenceSeq::getSupplementSubRef(std::string& sub_ref, std::string& ssuplmt)
{
	supplementRef(sub_ref, ssuplmt);
}

void ReferenceSeq::getSupplementSubRef(std::string chrom,int pos, int length, std::string& ssuplmt)//return supplement sub-reference.
{
	std::string temp=this->fp.parseFasta(chrom,pos,length);
	supplementRef(temp, ssuplmt);
}

void ReferenceSeq::setPath(std::string path)
{
	this->path=path;
	fp.setPath(path);
}

void ReferenceSeq::mapChromIDName()//map chrom ID with Name, by using XX.fasta.fai file
{
	fp.mapChromIDName();
}

void ReferenceSeq::loadChromIDName()
{
	fp.loadChromIDName();
}

void ReferenceSeq::getChromIDByName(std::string name, int& id)//given chrom name, return chrom id.
{
	int size=fp.vid_name.size();
	for(int i=0;i<size;i++)
	{
		if(fp.vid_name[i].second==name)
			id=fp.vid_name[i].first;
	}
}

void ReferenceSeq::getChromNameByID(int id, std::string& name)//given chrom id, return chrom name.
{
	int size=fp.vid_name.size();
	
	int begin=0;
	int end=size-1;
	while(begin<end)
	{
		int mid=(end-begin)/2;
		if(fp.vid_name[mid].first==id)
		{
			name=fp.vid_name[mid].second;
			return;
		}
		else if(fp.vid_name[mid].first<id)
		{
			begin=mid;
		}
		else
		{
			end=mid;
		}

		if(end-begin==1)
		{
			if(fp.vid_name[begin].first==id)
			{
				name=fp.vid_name[begin].second;
				return;
			}
			else
			{
				name=fp.vid_name[end].second;
				return;
			}
		}
	}//end of while
}

//-----------------private functions--------------------------------------------------------

/*
Descrition:
	Convert a sub ref to it's supplementary ref.
Input:
	sorgnl: string, a sub-ref;
Output:
	ssumplmt: string, a supplementary sub-ref. 
*/
void ReferenceSeq::supplementRef(std::string& sorgnl, std::string& ssuplmt)//convert to summplement ref
{	
	int len=sorgnl.length();
	ssuplmt="";
	for(int i=0;i<len;i++)
	{
		if(sorgnl[i]=='A' || sorgnl[i]=='a')
			ssuplmt+="T";
		else if(sorgnl[i]=='C' || sorgnl[i]=='c')
			ssuplmt+="G";
		else if(sorgnl[i]=='G' || sorgnl[i]=='g')
			ssuplmt+="C";
		else if(sorgnl[i]=='T' || sorgnl[i]=='t')
			ssuplmt+="A";
		else
			ssuplmt+=sorgnl[i];
	}
}

//convert LowerCase to Upper 
void ReferenceSeq::cvt2Uppercase()
{
	
}
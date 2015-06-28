#include<string>
#include<vector>
#include<algorithm>
#include<fstream>
#include<iostream>
#include"TrainingSet.h"
#include"public_parameters.h"

using namespace std;

void TrainingSet::setFilePath(std::string fbcmk, std::string ffeature, std::string ftrain)
{
	this->fbcmk=fbcmk;
	this->ffeature=ffeature;
	this->ftrain=ftrain;
}

void TrainingSet::gnrtTrainingSet() 
{
	prepareBcmk();//prepare training data

	int pos;
	double clip;
	double sclip_map;
	double discor;
	double cov_diff;
	string dir;

	ofstream fout_train;
	fout_train.open(ftrain.c_str());

	//read in feature data 
	ifstream fin_feature;
	fin_feature.open(ffeature.c_str());
	int vsize=this->vbcmk.size();
	
	//memset(vbcmk_exist,0, vsize);
	vector<int> vcandidates;
	vcandidates.clear();
	int bexist;
	while(fin_feature>>pos>>clip>>cov_diff>>sclip_map>>discor>>dir)
	{
		bexist=0;
		//check whether exist in benchmark
		for(int i=0;i<vsize;i++)
		{
			int lrange=pos-BASE_SLACK;
			int rrange=pos+BASE_SLACK;
			if(vbcmk[i]>=lrange && vbcmk[i]<=rrange)
			{
				bexist=1;
				vcandidates.push_back(pos);
				break;
			}
		}
		fout_train<<bexist<<" 1:"<<clip<<" 2:"<<cov_diff<<" 3:"<<sclip_map<<" 4:"<<discor<<endl;
	}
	
	//output false negative 
	ofstream fout_fn;
	fout_fn.open("false_negative.txt");
	fout_fn<<"test"<<endl;
	int vcan_size=vcandidates.size();
	for(int i=0;i<vsize;i++)
	{
		bexist=0;
		for(int j=0;j<vcan_size;j++)
		{
			int pos=vcandidates[j];
			int lrange=pos-BASE_SLACK;
			int rrange=pos+BASE_SLACK;
			if(vbcmk[i]>=lrange && vbcmk[i]<=rrange)
			{
				bexist=1;
				break;
			}
		}
		if(bexist==0)
			fout_fn<<vbcmk[i]<<endl;
	}

	fout_fn.close();
	fin_feature.close();
	fout_train.close();
}


void TrainingSet::gnrtTestingSet()
{
	int pos;
	double clip;
	double sclip_map;
	double discor;
	double cov_diff;
	string dir;

	ofstream fout_test;
	fout_test.open(ftrain.c_str());

	ifstream fin_feature;
	fin_feature.open(ffeature.c_str());
	while(fin_feature>>pos>>clip>>cov_diff>>sclip_map>>discor>>dir)
	{
		fout_test<<"1:"<<clip<<" 2:"<<cov_diff<<" 3:"<<sclip_map<<" 4:"<<discor<<endl;
	}

	fout_test.close();
	fin_feature.close();
}

bool TrainingSet::binarySearch(int pos)
{
	return false;
}

void TrainingSet::prepareBcmk()
{
	ifstream fin;
	fin.open(fbcmk.c_str());
	
	string chrom;
	int start, end;
	vector<int> vtemp;
	vtemp.clear();
	while(fin>>chrom>>start>>end) //benchmark format: chrom_name, start, end.
	{
		vtemp.push_back(start);
		vtemp.push_back(end);
	}
	std::sort(vbcmk.begin(),vbcmk.end()); //sort according positions

	//remove repeat ones 
	int vsize=vtemp.size();
	int last=vtemp[0];
	int cur;
	for(int i=1;i<vsize;i++)
	{
		cur=vtemp[i];
		if(cur==last)
		{
			continue;
		}
		else
		{
			this->vbcmk.push_back(last);
			last = cur;
		}
	}
	this->vbcmk.push_back(last);
	fin.close();
}
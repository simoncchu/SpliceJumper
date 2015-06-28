#include<iostream>
#include<cstdlib>
#include<fstream>
#include"bam_parse.h"
#include"fasta_parser.h"
//#include"ResultAnalyze/result_analyzation.h"
#include"fai_parser.h"
#include"algorithms/local_alignment.h"
#include"Coverage.h"
#include"CandidateSitesCaller.h"
#include"TrainingSet.h"
#include"HardClipReads.h"

void Usage()
{
	std::cout<<"Usage:"<<std::endl;
	std::cout<<"Rmjf file1 file2  ==> default cutoff=3"<<std::endl;
	std::cout<<"Rmjf -c value file1 file2  ==> set cutoff=value"<<std::endl;
}


int main(int argc, char* argv[])
{
	BamParse bp(argv[1]);
	bp.loadIndex();
	int refsize=bp.getRefInfo();
	
	//Run chrom by chrom
	//for(int cid=0;cid<refsize;cid++)
	for(int cid=0;cid<1;cid++)
	{
		cid=10; //chr1
		//cid=1;//chr11
		//find out all the candidate sitess 

		string chrom_name;
		int chrom_len;
		bp.getChromNameLength(cid, chrom_name, chrom_len);
		double mean_insert=250;
		double devi_insert=75;

		cout<<chrom_name<<" "<<chrom_len<<endl;//////////////////////////////////////////////////////////////////////////////////////////////////////
		Coverage cov; //step1: calc coverage.
		cov.initFiles();
		cov.setInsertSize(mean_insert,devi_insert);
		cov.pileup(cid,bp);

		//add original raw hard-clip read to clip_reads.txt////////////////////////////////////////////////////////////////////////////////////////
		HardClipReads hcr;
		hcr.traceRawReadById(argv[6],argv[7]);
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		CandidateSitesCaller csc;
		csc.setChromNameLen(chrom_name, chrom_len);
		csc.setInsertSize(mean_insert,devi_insert);
		csc.setRefPath((string)argv[2]);
		csc.callCandidateSites();
		csc.outputFeatures((string)argv[3]);
		csc.outputGraph();
		
		string sbrkpnt=(string)argv[3];
		sbrkpnt+="_brkpnt.txt";
		TrainingSet ts;
		ts.setFilePath((string)argv[4],sbrkpnt,(string)argv[5]);
		ts.gnrtTrainingSet();
		//ts.gnrtTestingSet();
		
	}
	
	return 0;
}

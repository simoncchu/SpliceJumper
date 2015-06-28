#include<string>

const int MAX_LEN_EACH_LINE_FAI=100;// the maximum length of each line in a .fai file.
const std::string CHROM_ID_NAME="ChromIDName.txt";//file to save chrom id and name 

// when partition ref, record the position of pileup file for each candidate exon into a file.
const std::string FILE_EXON_POSCOV_PILEUP="exon_pos_in_pileup_file.txt";

const int FRONT_DIRECTION=1;
const int TAIL_DIRECTION=2;
const std::string CANDIDATE_EXON_SITES="candiate_exon_all.site";//file to save called candidate exon sites 

const std::string UNMAPPED_NO_EXON="unmapped_reads_no_exon.sam";//"ReadsCollection.cpp" save all the unmapped reads no exon contain
const std::string CLIP_NO_REGION="./region_reads/clipped_reads_no_region.txt";

const int EXTEND_EXON=100;// for each candidate exon, extend the range for 100 each both end. "ReAlignment.cpp"
const int EXTEND_REGION_HIGH=35;
const int EXTEND_REGION_LOW=20;
//phase1 candidate exon file name style: chromName_exonfilename PartitionRef.h 

const int READ_LENGTH=100;

const int READ_TYPE_FULLMAP=1;
const int READ_TYPE_CLIP=2;
const int READ_TYPE_LEFT_SOFTCLIP=21;
const int READ_TYPE_RIGHT_SOFTCLIP=22;
const int READ_TYPE_BOTH_SOFTCLIP=23;
const int READ_TYPE_LEFT_HARDCLIP=24;
const int READ_TYPE_RIGHT_HARDCLIP=25;
const int READ_TYPE_BOTH_HARDCLIP=26;
const int READ_TYPE_UNMAP=3;
const int READ_TYPE_OTHER=4;

const int READ_PAIR_MAP_TYPE_11=11;//both mapped 
const int READ_PAIR_MAP_TYPE_10=10;//read map, mate unmap
const int READ_PAIR_MAP_TYPE_01=12;//read unmap, mate map
const int READ_PAIR_MAP_TYPE_00=0;//both unmapped 

const double COV_CUTOFF=0.65;//(A-B)^2/A^2, normalized 

const int CUTOFF_LA_SGMT=7;//local alignment segment cutoff value, smaller than this value, then will not consider
const double CUTOFF_LA_OPTIMAL_PERCENT=0.8;// 80% percent of a segment is local aligned, then see this as a good match 

const int SMALLEST_INTRON=70;

const int CUTOFF_CALLINTRON_BY_READ=5; //if more than 5 reads support, then call this a INTRON 

const int NEIGHBOR_REGIONS=5;//JunctionCaller.cpp


//-----------------Added 07/21/14-----------------------------------------------------------------------------------
const std::string fname_cov="coverage.txt";
const std::string fname_fm="full_map_reads.txt";
const std::string fname_discor="discordant_reads.txt";//save all the discordant reads
const std::string fname_clip="clip_reads.txt";
const std::string fname_um="unmap_reads.txt";
const std::string fname_others="other_reads.txt";
const std::string fname_brkpnts="candidate_brkpnts.txt";
const std::string fname_graph="brkpnt_graph.txt";
const int MIN_CLIP_READS=1;//minimum supported reads to call out a candidate site
const int NO_CLIP_POS=1000;//indicates no clip happen
const int MAX_INTRON_SIZE=1500000;
const int BASE_SLACK=15;
const double HIT_LOCAL_ALIGNMENT=0.8;//ratio to consider a local alignment is a hit 
const int CANON_SLACK=3;

const int LDISCORDANT=0;
const int RDISCORDANT=1;
const int CONCORDANT=2;
const int WRONG_PAIR=3;

//tmp files 
const std::string fname_icp="is_clip_pos.txt";
const std::string fname_ldiscor="ldiscor.txt";
const std::string fname_rdiscor="rdiscor.txt";
const std::string fname_hard_clip="hard_clip.txt";
const std::string fname_left_hclip="left_hard_clip_id.txt";
const std::string fname_right_hclip="right_hard_clip_id.txt";
const std::string fname_lhclip_raw="left_hard_clip_raw_reads.txt";
const std::string fname_rhclip_raw="right_hard_clip_raw_reads.txt";
BAMTOOLS=../bamtools-master/include
BAMTOOLS_LD=../bamtools-master/lib

CC = g++
CFLAGS =  -O3 -Wall -I$(BAMTOOLS) -L$(BAMTOOLS_LD) -Wl,-rpath,$(BAMTOOLS_LD)

all: SpliceJumper

public_func.o :	public_func.cpp public_func.h
	$(CC) $(CFLAGS) -c public_func.cpp

Alignment.o: Alignment.cpp Alignment.h
	$(CC) $(CFLAGS) -c Alignment.cpp
	
bam_parse.o : bam_parse.cpp bam_parse.h
	$(CC) $(CFLAGS) -c bam_parse.cpp

Coverage.o: Coverage.cpp Coverage.h
	$(CC) $(CFLAGS) -c Coverage.cpp

fai_parser.o: fai_parser.h
	$(CC) $(CFLAGS) -c fai_parser.cpp

fasta_parser.o: fasta_parser.cpp fasta_parser.h
	$(CC) $(CFLAGS) -c fasta_parser.cpp

ReferenceSeq.o : ReferenceSeq.cpp ReferenceSeq.h
	$(CC) $(CFLAGS) -c ReferenceSeq.cpp

local_alignment.o: algorithms/local_alignment.cpp algorithms/local_alignment.h
	$(CC) $(CFLAGS) -c algorithms/local_alignment.cpp

CandidateSitesCaller.o: CandidateSitesCaller.cpp CandidateSitesCaller.h
	$(CC) $(CFLAGS) -c CandidateSitesCaller.cpp

TrainingSet.o: TrainingSet.cpp TrainingSet.h
	$(CC) $(CFLAGS) -c TrainingSet.cpp
	
HardClipReads.o: HardClipReads.cpp HardClipReads.h
	$(CC) $(CFLAGS) -c HardClipReads.cpp

main.o: main.cpp 
	$(CC) $(CFLAGS) -c main.cpp

SpliceJumper: public_func.o Alignment.o bam_parse.o Coverage.o fai_parser.o fasta_parser.o ReferenceSeq.o local_alignment.o CandidateSitesCaller.o TrainingSet.o HardClipReads.o main.o 
	$(CC) $(CFLAGS) -o JunctionClassifier_1 public_func.o Alignment.o bam_parse.o Coverage.o fai_parser.o fasta_parser.o ReferenceSeq.o local_alignment.o CandidateSitesCaller.o TrainingSet.o HardClipReads.o main.o \
	-lbamtools -lz -lm
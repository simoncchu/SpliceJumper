
Usage: ./SpliceJumper [options] 

Required options:
		  -P/G        indicates collecting features / indicates output results 
                  -r FILE   reference file(indexed)
                  -b/i FILE   input bam/prediction file(index and sorted)
                  -o FILE   output file name
                  -l INT    read length
Optimal options
                  -c INT    slack value for split position with default 3
                  -m DOUBLE mean insert size
                  -v DOUBLE standard variation of insert size
                  
1.Feature collection:
./SpliceJumper -P -r ./human_g1k_v37.fasta -b ./simulated_data_input.bam -o test_sim.txt -l 100 -c 5

2. Training and prediction
2.1 Training
python easy.py training_data

2.2 Predicting
..\windows\svm-predict data.scale Trained.model predictResult

3. Output generated splicing junctions
./SpliceJumper -G -i prediction.txt -o result.txt
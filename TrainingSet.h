#ifndef _H_TESGNINIART_
#define _H_TESGNINIART_

#include<string>
#include<vector>

class TrainingSet
{
public:
	TrainingSet(){vbcmk.clear();}
	void setFilePath(std::string fbcmk, std::string ffeature, std::string ftrain);

public:
	void gnrtTrainingSet();
	void gnrtTestingSet();

private:
	void prepareBcmk();
	bool binarySearch(int pos);
private:
	std::string ftrain;
	std::string fbcmk;
	std::string ffeature;
	std::vector<int> vbcmk;
};

#endif


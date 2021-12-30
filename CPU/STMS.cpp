#include <iostream>
#include <stdint.h>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <algorithm>
#include "MurmurHash3.h"
using namespace std;

#define WORD_TYPE uint64_t //actual type of a word
#define KEY_SIZE 4//bytes num of flowId/eleId

class STM_DF
{
private:
	static const unsigned int wordLen = sizeof(WORD_TYPE) * 8;
	unsigned int k;
	unsigned int myWordNum;
	WORD_TYPE* myWordArray;
public:
	STM_DF() {
		k = 0;
		myWordNum = 0;
		myWordArray = NULL;
	}
	STM_DF(unsigned int wordNum, unsigned int hashNum) {
		if (hashNum > 8) {
			cout << "STM_DF: sorry, in this implementation, we can not support more than 8 hash funcitons, you need to modifiy the insert function!";
			exit(-1);
		}
		myWordNum = wordNum;
		myWordArray = new WORD_TYPE[myWordNum];
		k = hashNum;

		int tempNum = ceil(((double)sizeof(WORD_TYPE)) / sizeof(unsigned int));

		//initialize
		for (unsigned int i = 0; i < myWordNum; i++) {
			myWordArray[i] = 0;
			for (int j = 0; j < tempNum; j++) {
				myWordArray[i] = myWordArray[i] << 32;
				myWordArray[i] += 0xAAAAAAAA;
			}
		}
	}
	~STM_DF() {
		if (myWordArray) {
			delete[]myWordArray;
		}
	}
	// return true: the element has been recorded before this insertion
	bool insert(uint64_t* hashValue) {

		//use the 32th bit as the status (count from the left)
		bool status = (hashValue[0] >> 32) & 1;

		//we use the second 32 bit for finding the word index
		WORD_TYPE* aWord = &(this->myWordArray[(hashValue[0] & 0xffff) % this->myWordNum]);

		//set the k bits
		//we use the remaining 64 bits to find the k bits in a word
		//we suppose a word will not longer than 255 bits, thus 64 bits can be divided into 8 8-bit values for at most 8 hash functions
		bool flag = true;
		WORD_TYPE mask = 0;
		unsigned int bitIndex = 0;
		for (unsigned int i = 0; i < this->k; i++) {
			bitIndex = (hashValue[1] >> (i << 3)) % this->wordLen;
			if (flag) {
				bool findedBit = ((*aWord) >> bitIndex) & 1;//the rightmost bit in a word is considered as word[0]
				if (findedBit != status) {
					flag = false;
				}
			}
			if (!flag) {
				if (status) {//set to 1
					*aWord = (*aWord) | (((WORD_TYPE)1) << bitIndex);
				}
				else {//set to 0
					*aWord = (*aWord) & (~((WORD_TYPE)1 << bitIndex));
				}
			}
		}
		return flag;
	}
};

//for recording information in the unordered_map and unordered_set
struct Cmp {
	bool operator()(const char* a, const char* b) const {
		return memcmp(a, b, KEY_SIZE) == 0;
	}
};
struct HashFunc {
	unsigned int operator()(const char* key) const {
		unsigned int hashValue = 0;
		MurmurHash3_x86_32(key, KEY_SIZE, 0, &hashValue);
		return hashValue;
	}
};

class STMS
{
private:
	STM_DF myDF;
	double p;
	unordered_map<char*, unordered_set<char*, HashFunc, Cmp>, HashFunc, Cmp> recordingInfo;//simulate off-chip recording
public:
	STMS(unsigned int bitsNum, unsigned int hashNum, double samplingRatio) :myDF(bitsNum / sizeof(WORD_TYPE), hashNum) {
		p = samplingRatio;
	}
	
	void sampleFlowElements(char* flowId, char* eleId) {
		char key[KEY_SIZE];
		for (int i = 0; i < KEY_SIZE; i++) {
			key[i] = flowId[i] ^ eleId[i];
		}

		uint64_t hashValue[2]{ 0 };
		MurmurHash3_x64_128(key, KEY_SIZE, 0, &hashValue);

		//since sampling ratio will not be too small, 31 bits is enough for pre-sampling
		if ((unsigned int)(hashValue[0] >> 33) < (this->p * INT32_MAX)) {
			bool hasRecordedFlag = this->myDF.insert(hashValue);//flag==true means STM_DF has recorded this key, and we should not send the element to off-chip

			if (!hasRecordedFlag) {
				this->recordingInfo[flowId].insert(eleId);
			}
		}
	}

	unsigned int estimate(char* flowId, double p_e) {
		if (this->recordingInfo.find(flowId) == recordingInfo.end()) {//flowId has not been recorded
			return 0;
		}
		else {
			unsigned int c_f = recordingInfo[flowId].size();
			unsigned int estimatedSpread = c_f / p_e;
			return estimatedSpread;
		}
	}
};

//save the result int a txt file
//each line is: realSpread estimatedSpread
void saveResults(string outputFilePath, STMS& stms, double pe, unordered_map<char*, unordered_set<char*, HashFunc, Cmp>, HashFunc, Cmp>& realFlowInfo) {
	ofstream fout;
	vector<unsigned int> realSpreads, estimatedSpreads;

	auto iter = realFlowInfo.begin();
	fout.open(outputFilePath, ios::out);
	unsigned int i = 0;
	while (fout.is_open() && iter != realFlowInfo.end()) {
		if (iter != realFlowInfo.begin()) {
			fout << endl;
		}
		unsigned int realSpread = (iter->second).size();
		unsigned int estimatedSpread = stms.estimate(iter->first, pe);

		realSpreads.push_back(realSpread);
		estimatedSpreads.push_back(estimatedSpread);

		fout << realSpread << " " << estimatedSpread;
		iter++;
	}
	if (!fout.is_open()) {
		cout << outputFilePath << " closed unexpectedlly";
	}
	else {
		fout.close();
	}
	sort(realSpreads.begin(), realSpreads.end(), greater<unsigned int>());
	cout << "";
}

void processPackets(STMS& stms, vector<pair<char*, char*>>& dataset) {
	//int count = 0;
	clock_t start = clock();
	for (unsigned int i = 0; i < dataset.size(); i++) {
		stms.sampleFlowElements(dataset[i].first, dataset[i].second);
	}

	clock_t current = clock();
	cout << dataset.size() << " lines: have used " << ((double)current - start) / 1000.0 << " seconds" << endl;

	double throughput = (dataset.size() / 1000000.0) / (((double)current - start) / 1000.0);
	cout << "throughput: " << throughput << "Mpps" << endl;
}

//calculate the probability p_e for estimating
double getPe(double p1, unsigned int k, unsigned int m, double tau, double c) {
	double fp = pow(2, -1 * (double)k);
	double fn_step1 = pow(1.0 - (double)k / m, p1 * tau);
	double fn_step2 = (1.0 + fn_step1) / 2.0;
	double fn = 1.0 - pow(fn_step2, k);
	double pe = p1 * (1 - fp * pow(1 - fn, c + 1));
	return pe;
}

//save the data in memory and count the real flow info
void getDataSet(string dataDir, unsigned int numOfMinutes, vector<pair<char*, char*>>& dataset, unordered_map<char*, unordered_set<char*, HashFunc, Cmp>, HashFunc, Cmp>& realFlowInfo) {
	char dataFileName[20];
	ifstream fin;
	char* flowId;
	char* eleId;
	clock_t start = clock();
	for (unsigned int i = 0; i < numOfMinutes; i++) {
		sprintf_s(dataFileName, "%02d.dat", i);
		string oneDataFilePath = dataDir + string(dataFileName);
		fin.open(oneDataFilePath, ios::in | ios::binary);

		while (fin.is_open() && fin.peek() != EOF) {
			flowId = new char[KEY_SIZE] {0};
			eleId = new char[KEY_SIZE] {0};

			fin.read(flowId, KEY_SIZE);
			fin.read(eleId, KEY_SIZE);

			dataset.push_back(make_pair(flowId, eleId));
			realFlowInfo[flowId].insert(eleId);//count the real flow info

			if (dataset.size() % 5000000 == 0) {//output someting to check the procedure
				clock_t current = clock();
				cout << "have added " << dataset.size() << " packets, have used " << ((double)current - start) / 1000.0 << " seconds." << endl;
			}
		}

		if (!fin.is_open()) {
			cout << "dataset file" << dataFileName << "closed unexpectedlly";
			exit(-1);
		}
		else {
			fin.close();
		}
	}

	clock_t current = clock();
	cout << "have added " << dataset.size() << " packets, have used " << ((double)current - start) / 1000.0 << " seconds" << endl;

	//count the totol number of flows and distinct elements
	auto iter = realFlowInfo.begin();
	unsigned int totalSpread = 0;
	while (iter != realFlowInfo.end()) {
		totalSpread += (iter->second).size();
		iter++;
	}
	cout << "there are " << realFlowInfo.size() << " flows and " << totalSpread << " distinct elements" << endl;
}

int main()
{
	//prepare the dataset
	cout << "prepare the dataset" << endl;
	string dataDir = R"(G:\networkTrafficData\2016\processed_binary\)";
	unsigned int numOfMinutes = 1;
	vector<pair<char*, char*>> dataset;
	unordered_map<char*, unordered_set<char*, HashFunc, Cmp>, HashFunc, Cmp> realFlowInfo;
	getDataSet(dataDir, numOfMinutes, dataset, realFlowInfo);

	//prepare the sampling
	cout << endl;
	cout << "prepare the sampling" << endl;
	unsigned int bitsNum = 10000;
	unsigned int hashNum = 4;
	double samplingRatio = 0.2;
	STMS stms(bitsNum, hashNum, samplingRatio);

	//prepare the estimate parameters
	double tau = 62200;//get from historical traffic data
	double c = 44;
	double pe = getPe(samplingRatio, hashNum, bitsNum, tau, c);


	//start sampling
	cout << endl;
	cout << "start sampling" << endl;
	processPackets(stms, dataset);
	cout << "Done" << endl;

	//save the result in files
	cout << endl;
	cout << "save the result in spreads.txt ..." << endl;
	string outputFilePath = "./spreads.txt";
	saveResults(outputFilePath, stms, pe, realFlowInfo);

	//release resources
	cout << endl;
	cout << "release resources..." << endl;
	for (int i = 0; i < dataset.size(); i++) {
		delete[](dataset[i].first);
		delete[](dataset[i].second);
	}

	return 0;
}
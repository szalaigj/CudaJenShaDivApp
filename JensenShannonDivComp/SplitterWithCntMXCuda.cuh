#ifndef SPLITTER_WITH_CNT_MX_CUH_
#define SPLITTER_WITH_CNT_MX_CUH_

#include "SplitterWithCntMX.hpp"

class SplitterWithCntMXCuda : public ISplitterWithCntMX
{
public:
	SplitterWithCntMXCuda(std::string& sequence, ICountMatrixHolder& countMatrixHolder) :
		sequence(sequence), countMatrixHolder(countMatrixHolder)
	{
	}

	ICountMatrixHolder& getCountMatrixHolder()
	{
		return countMatrixHolder;
	}

	std::vector<std::string> split(int startIdx, int endIdx);
	double * computeDivergenceFromCountArray(BYTE * countArray, int startIdx, int endIdx, int seqLen);
	void determineMaxDivAndPos(int startIdx, int seqLen, double * divs_A, double * divs_C, double * divs_G, double * divs_T, double& maxDivergence, int& maxPos);
	double computeSignificance(int, double);
	std::vector<std::string> checkSignificance(int, int, int, double);

private:
	std::string& sequence;
	ICountMatrixHolder& countMatrixHolder;
};

#endif /* SPLITTER_WITH_CNT_MX_CUH_ */
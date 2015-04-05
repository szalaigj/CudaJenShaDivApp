// standard library includes
#include <iostream>
#include <math.h>

#include <gsl/gsl_cdf.h>

#include "SplitterWithCntMX.hpp"

extern const float betaParam;
extern const float aParam;
extern const float bParam;
extern const double significanceThreshold;
extern const int minSeqLength;
extern const std::string symbols;

std::vector<std::string> SplitterWithCntMX::split(int startIdx, int endIdx)
{
	double maxDivergence = -1.0;
	int seqLen = endIdx - startIdx + 1;
	int maxPos = -1;
	if (seqLen >= minSeqLength)
	{
		for (int pos = startIdx + 1; pos <= endIdx; pos++)
		{
			computeDivergenceForPos(startIdx, endIdx, pos, maxDivergence, maxPos);
		}
	}
	double significance = computeSignificance(endIdx - startIdx + 1, maxDivergence);
	std::vector<std::string> result = checkSignificance(startIdx, endIdx, maxPos, significance);
	return result;
}

void SplitterWithCntMX::computeDivergenceForPos(int startIdx, int endIdx, int pos, double& maxDivergence, int& maxPos)
{
	double a = countMatrixHolder.queryFrequency('A', startIdx, pos - 1);
	double c = countMatrixHolder.queryFrequency('C', startIdx, pos - 1);
	double g = countMatrixHolder.queryFrequency('G', startIdx, pos - 1);
	double t = countMatrixHolder.queryFrequency('T', startIdx, pos - 1);
	double frequenciesPrefix[4] = { a, c, g, t };
	a = countMatrixHolder.queryFrequency('A', pos, endIdx);
	c = countMatrixHolder.queryFrequency('C', pos, endIdx);
	g = countMatrixHolder.queryFrequency('G', pos, endIdx);
	t = countMatrixHolder.queryFrequency('T', pos, endIdx);
	double frequenciesPostfix[4] = { a, c, g, t };
	double weightPrefix = (double)(pos - startIdx) / (double)(endIdx - startIdx + 1);
	double weightPostfix = (double)(endIdx - pos + 1) / (double)(endIdx - startIdx + 1);
	double divergence = jenShaDivComputer.computeDivergence(frequenciesPrefix, frequenciesPostfix,
		weightPrefix, weightPostfix, symbols.length(), symbols.length());
	if (maxDivergence < divergence)
	{
		maxDivergence = divergence;
		maxPos = pos;
	}
}

double SplitterWithCntMX::computeSignificance(int N, double maxDivergence)
{
	double significance = 0.0;
	if (maxDivergence > 0.0)
	{
		double NEff = aParam * log(N) + bParam;
		significance = gsl_cdf_chisq_P(2 * N * log(2) * betaParam * maxDivergence, symbols.length() - 1);
		significance = pow(significance, NEff);
	}
	return significance;
}

std::vector<std::string> SplitterWithCntMX::checkSignificance(int startIdx, int endIdx, int maxPos, double significance)
{
	std::vector<std::string> result;
	if ((significance > significanceThreshold) && (maxPos >= 0))
	{
		std::cout << "New splitpoint on pos " << maxPos << std::endl;
		std::vector<std::string>& resultPrefix = split(startIdx, maxPos - 1);
		std::vector<std::string>& resultPostfix = split(maxPos, endIdx);
		result.resize(resultPrefix.size() + resultPostfix.size());
		result = std::vector<std::string>(resultPrefix);
		result.insert(result.end(), resultPostfix.begin(), resultPostfix.end());
	}
	else
	{
		result.resize(1);
		std::string resultSeq = sequence.substr(startIdx, endIdx - startIdx + 1);
		result.push_back(resultSeq);
	}
	return result;
}
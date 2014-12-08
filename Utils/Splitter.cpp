#include <math.h>

#include "Splitter.hpp"

extern const float betaParam;
extern const float aParam;
extern const float bParam;
extern const double significanceThreshold;
extern const int minSeqLength;
extern const std::string symbols;

std::vector<std::string>& Splitter::split(std::string& sequence)
{
	std::vector<std::string> result;
	double maxDivergence = -1.0;
	std::string seqPrefixForMax = "";
	std::string seqPostfixForMax = "";
	std::map <char, long> chrCountsPrefix;
	getFrequencyComputer().initChrCounts(chrCountsPrefix);
	std::map <char, long> chrCountsPostfix;
	getFrequencyComputer().initChrCounts(chrCountsPostfix);
	getFrequencyComputer().countCharsInSeq(sequence, chrCountsPostfix);
	int seqLen = sequence.length();
	if (seqLen >= minSeqLength)
	{
		for (int pos = 1; pos < seqLen; pos++)
		{
			computeDivergenceForPos(sequence, maxDivergence, seqPrefixForMax, seqPostfixForMax,
				pos, chrCountsPrefix, chrCountsPostfix);
		}
	}
	double significance = computeSignificance(sequence, maxDivergence);
	result = checkSignificance(sequence, seqPrefixForMax, seqPostfixForMax, significance);
	return result;
}

void Splitter::computeDivergenceForPos(std::string& sequence, double& maxDivergence,
	std::string& seqPrefixForMax, std::string& seqPostfixForMax, int pos,
	std::map <char, long>& chrCountsPrefix, std::map <char, long>& chrCountsPostfix)
{
	int seqLength = sequence.length();
	std::string sequencePrefix = sequence.substr(0, pos);
	std::string sequencePostfix = sequence.substr(pos);
	char currentChar = sequence[pos - 1];
	frequencyComputer.increaseCountChar(currentChar, chrCountsPrefix);
	frequencyComputer.decreaseCountChar(currentChar, chrCountsPostfix);
	double * frequenciesPrefix = frequencyComputer.determineFrequencies(chrCountsPrefix, seqLength);
	double * frequenciesPostfix = frequencyComputer.determineFrequencies(chrCountsPostfix, seqLength);
	double weightPrefix = (double)pos / (double)seqLength;
	double weightPostfix = ((double)seqLength - (double)pos) / (double)seqLength;
	double divergence = jenShaDivComputer.computeDivergence(frequenciesPrefix, frequenciesPostfix,
		weightPrefix, weightPostfix, symbols.length(), symbols.length());
	if (maxDivergence < divergence)
	{
		maxDivergence = divergence;
		seqPrefixForMax = sequencePrefix;
		seqPostfixForMax = sequencePostfix;
	}
}

double Splitter::computeSignificance(std::string& sequence, double maxDivergence)
{
	double significance = 0.0;
	if (maxDivergence > 0.0)
	{
		int N = sequence.length();
		double NEff = aParam * log(N) + bParam;
		//chi_squared chiSquared(symbols.length() - 1);
		//significance = boost::math::cdf(chiSquared, N * log(2) * betaParam * maxDivergence);
		significance = pow(significance, NEff);
	}
	return significance;
}

std::vector<std::string>& Splitter::checkSignificance(std::string& sequence, std::string& seqPrefixForMax,
	std::string& seqPostfixForMax, double significance)
{
	std::vector<std::string> result;
	if (significance > significanceThreshold)
	{
		std::vector<std::string>& resultPrefix = split(seqPrefixForMax);
		std::vector<std::string>& resultPostfix = split(seqPostfixForMax);
		result.resize(resultPrefix.size() + resultPostfix.size());
		result = std::vector<std::string>(resultPrefix);
		result.insert(result.end(), resultPostfix.begin(), resultPostfix.end());
	}
	else
	{
		result.resize(1);
		result.push_back(sequence);
	}
	return result;
}
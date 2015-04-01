#include "JenShaDivComputer.hpp"

double JenShaDivComputer::computeDivergence(double * frequencies1, double * frequencies2, double weight1, double weight2, int freq1Length, int freq2Length)
{
	if (freq1Length != freq2Length)
	{
		throwErrorMsg();
	}
	double * weightedFrequencies1 = computeWeightedFrequencies(frequencies1, weight1, freq1Length);
	double * weightedFrequencies2 = computeWeightedFrequencies(frequencies2, weight2, freq2Length);
	double * sumOfFrequencies = computeSum(weightedFrequencies1, weightedFrequencies2, freq1Length);
	double result = getEntropyComputer().computeEntropy(sumOfFrequencies, freq1Length)
		- weight1 * getEntropyComputer().computeEntropy(frequencies1, freq1Length)
		- weight2 * getEntropyComputer().computeEntropy(frequencies2, freq2Length);
	// The followings are needed otherwise there is memory leak.
	delete weightedFrequencies1;
	delete weightedFrequencies2;
	delete sumOfFrequencies;
	return result;
}
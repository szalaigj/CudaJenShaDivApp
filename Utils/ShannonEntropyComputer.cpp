#include "ShannonEntropyComputer.hpp"

double ShannonEntropyComputer::computeEntropy(double * frequencies, int size)
{
	double entropy = 0.0;
	double frequency;
	for (int idx = 0; idx < size; idx++)
		{
		frequency = frequencies[idx];
		if (frequency != 0.0)
		{
			entropy -= frequency * log(frequency);
		}
	}
	return entropy;
}
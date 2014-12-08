#ifndef JEN_SHA_DIV_COMP_HPP_
#define JEN_SHA_DIV_COMP_HPP_

#include <stdexcept>

#include "ShannonEntropyComputer.hpp"

class IJenShaDivComputer
{
	public:
		virtual double computeDivergence(double * frequencies1, double * frequencies2, double weight1, double weight2, int freq1Length, int freq2Length) = 0;
};

class JenShaDivComputer : public IJenShaDivComputer
{
	public:
		JenShaDivComputer(ShannonEntropyComputer& entropyComputer) : entropyComputer(entropyComputer)
		{
		}

		ShannonEntropyComputer& getEntropyComputer()
		{
			return entropyComputer;
		}

		virtual double computeDivergence(double *, double *, double, double, int, int);

	private:
		void throwErrorMsg()
		{
			std::string errmsg("The two frequency inputs have not same length.");
			throw std::runtime_error(errmsg);
		}

		double * computeWeightedFrequencies(double * frequencies, double weight, int freqLength)
		{
			double * weightedFrequencies = new double[freqLength];
			for (int idx = 0; idx < freqLength; idx++)
			{
				weightedFrequencies[idx] = frequencies[idx] * weight;
			}
			return weightedFrequencies;
		}

		double * computeSum(double * frequencies1, double * frequencies2, int freqLength)
		{
			double * sumOfFrequencies = new double[freqLength];
			for (int idx = 0; idx < freqLength; idx++)
			{
				sumOfFrequencies[idx] = frequencies1[idx] + frequencies2[idx];
			}
			return sumOfFrequencies;
		}

		ShannonEntropyComputer& entropyComputer;
};

#endif /* JEN_SHA_DIV_COMP_HPP_ */
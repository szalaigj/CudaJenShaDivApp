#ifndef SHA_ENT_COMP_HPP_
#define SHA_ENT_COMP_HPP_

// standard includes
#include <math.h>

class IShannonEntropyComputer
{
	public:
		virtual double computeEntropy(double * frequencies, int size) = 0;
};

class ShannonEntropyComputer : public IShannonEntropyComputer
{
	public:
		double computeEntropy(double *, int);
};

#endif /* SHA_ENT_COMP_HPP_ */
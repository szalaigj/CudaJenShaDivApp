#ifndef CHISQUARE_HPP_
#define CHISQUARE_HPP_

#include "UIncGamma.cpp"

class IChiSquaredCDFComputer
{
	public:
		virtual double computeValue(int k, double x) = 0;
};

class ChiSquaredCDFComputer : public IChiSquaredCDFComputer
{
	public:
		ChiSquaredCDFComputer(UIncGamma& uIncGamma) : uIncGamma(uIncGamma)
		{
		}

		double computeValue(int k, double x);
    private:
		UIncGamma& uIncGamma;
};
#endif /* CHISQUARE_HPP_ */
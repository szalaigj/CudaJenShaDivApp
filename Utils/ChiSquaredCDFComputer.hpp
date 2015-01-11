#ifndef CHISQUARE_HPP_
#define CHISQUARE_HPP_

class IChiSquaredCDFComputer
{
	public:
		virtual double gammaApprox(double s) = 0;
		virtual double lowerIncompleteGamma(double s, double x) = 0;
		virtual double chiSquared(int k, double x) = 0;
};

class ChiSquaredCDFComputer : public IChiSquaredCDFComputer
{
	public:
		double gammaApprox(double s);
		double lowerIncompleteGamma(double s, double x);
		double chiSquared(int k, double x);
};
#endif /* CHISQUARE_HPP_ */
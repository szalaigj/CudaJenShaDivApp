#include <math.h>

#include "ChiSquaredCDFComputer.hpp"

# define PI 3.14159265358979323846264338

// The following implementation uses Gergo Nemes's proposition for approximation:
//  Nemes, Gergõ (2010)
//  "New asymptotic expansion for the Gamma function"
//  Archiv der Mathematik 95 (2): 161–169
double ChiSquaredCDFComputer::gammaApprox(double s)
{
	double temp = (1.0 / exp(1.0)) * (s + 1.0/(12*s - 1.0/(10*s)));
	return sqrt(2 * PI / s) * pow(temp, s);
}

// The following is not the best implementation
// because for s = 1.5 and x = 17976.7740590471075
// the result is -1.#IND000000000000 instead of 0.88622692545275801365
double ChiSquaredCDFComputer::lowerIncompleteGamma(double s, double x)
{
	int sumRange = 100;
	double sum = 0.0;
	for (int idx = 0; idx <= sumRange; idx++)
	{
		double temp = pow(x, idx);
		double prod = 1.0;
		for (int subIdx = 0; subIdx <= idx; subIdx++)
		{
			prod *= (s + subIdx);
		}
		temp = temp / prod;
		sum += temp;
	}
	return exp(-x) * pow(x, s) * sum;
}

double ChiSquaredCDFComputer::chiSquared(int k, double x)
{
	double s = (double)k / 2.0;
	double halfX = x / 2.0;
	//double result = lowerIncompleteGamma(s, halfX) / gammaApprox(s);
	int sumRange = 100;
	double sum = 0.0;
	for (int idx = s; idx <= sumRange; idx++)
	{
		double y = log(halfX);
		y *= idx;
		y = y - log(gammaApprox(idx + 1));
		sum += exp(y);
	}
	double result = exp(-halfX) * sum;
	return result;
}
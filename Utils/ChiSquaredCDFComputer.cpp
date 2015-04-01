#include <math.h>

#include "ChiSquaredCDFComputer.hpp"

double ChiSquaredCDFComputer::computeValue(int k, double x)
{
	double s = (double)k / 2.0;
	double halfX = x / 2.0;
	double result = 1 - exp(log(uIncGamma.computeValue(s, halfX)) - log(tgamma(s)));
	return result;
}
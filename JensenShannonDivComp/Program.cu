/*
* This implementation is based on the article:
*	Grosse, I., Bernaola-Galvan, P., Carpena, P., Roman-Roldan, R., Oliver, J., & Stanley, H. E. (2002).
*	Analysis of symbolic sequences using the Jensen-Shannon divergence.
*	Physical Review E, 65(4), 041905
*/

// standard library includes
#include <iostream>

// local CudaJenShaDivApp includes
#include "CudaBasicIncludes.cuh"
#include "../Utils/ShannonEntropyComputer.cpp"
#include "../Utils/FrequencyComputer.cpp"
#include "../Utils/JenShaDivComputer.cpp"
#include "../DataHandlingUtil/LoadData.cpp"
#include "../DataHandlingUtil/SubSeqWriter.cpp"
#include "../Utils/Splitter.cpp"
#include "../Utils/ChiSquaredCDFComputer.cpp"

//#include <gsl/gsl_sf_gamma.h>
//#include <gsl/gsl_randist.h>
//#include <gsl/gsl_cdf.h>

// used symbols in sequences
const std::string symbols = "ACGT";

// the following values are referred to in article (Table 1)
// where the number of symbols is 4
const float betaParam = 0.8f;
const float aParam = 2.44f;
const float bParam = -6.15f;
const float hUsedParams[3] = { betaParam, aParam, bParam };

const double significanceThreshold = 0.9;
const int minSeqLength = 13;

__constant__ int sequenceLength;
__constant__ double usedParams[3];

int main(int argc, char** argv)
{
	BaseDataLoader loader;
	ShannonEntropyComputer entropyComputer;
	FrequencyComputer frequencyComputer(symbols);
	JenShaDivComputer jenShaDivComputer(entropyComputer);
	ChiSquaredCDFComputer chiSquaredCDFComputer;
	Splitter splitter(frequencyComputer, jenShaDivComputer, chiSquaredCDFComputer);
	BaseSubSeqWriter subSeqWriter;
	std::string filename(argv[1]);
	std::string inputData = loader.loadData(filename);
	int len = inputData.length();
	//std::cout << inputData << std::endl;
	cudaCheckError(cudaMemcpyToSymbol(sequenceLength, &len, sizeof(int)));
	cudaCheckError(cudaMemcpyToSymbol(usedParams, hUsedParams, 3 * sizeof(double)));
	//std::string seq = "ACAGCAGGATATGTCGTGCT";
	std::string seq = "ACACACACACACACGTGTGTGTGTGTG";
	double * frequencies = frequencyComputer.computeFrequency(&seq);
	double entropy = entropyComputer.computeEntropy(frequencies, symbols.length());
	std::cout << entropy << std::endl;
	std::string sequencePrefix = seq.substr(0, 5);
	std::string sequencePostfix = seq.substr(5);
	double * frequenciesPrefix = frequencyComputer.computeFrequency(&sequencePrefix);
	double * frequenciesPostfix = frequencyComputer.computeFrequency(&sequencePostfix);
	double weightPrefix = (double)sequencePrefix.length() / (double)seq.length();
	double weightPostfix = (double)sequencePostfix.length() / (double)seq.length();
	double divergence = jenShaDivComputer.computeDivergence(frequenciesPrefix, frequenciesPostfix,
		weightPrefix, weightPostfix, symbols.length(), symbols.length());
	std::cout << divergence << std::endl;
	//double d3 = gsl_cdf_chisq_P(35953.548118094215, 3.0);
	//double d2 = log(gsl_sf_gamma_inc_P(1.5, 17976.774059047108)) + gsl_sf_lngamma(17976.774059047108);
	//double d1 = gsl_sf_gamma(17976.774059047108) - gsl_sf_gamma_inc(1.5, 17976.774059047108);
	//double x = 2 * 547496 * log(2) * betaParam * 0.059212766850595649;
	//double temp = (double)(symbols.length() - 1) / 2.0;
	//double res1 = chiSquaredCDFComputer.lowerIncompleteGamma(temp, x / 2.0);
	//double res2 = chiSquaredCDFComputer.gammaApprox(temp);
	//double res = chiSquaredCDFComputer.chiSquared(symbols.length() - 1, 2 * 547496 * log(2) * betaParam * 0.059212766850595649);
	std::vector<std::string> subsequences = splitter.split(inputData);
	std::string outputFilename(argv[2]);
	subSeqWriter.writeData(outputFilename, &subsequences);
	std::cout << "Press any character and press enter to continue..." << std::endl;
	char chr;
	std::cin >> chr;
}
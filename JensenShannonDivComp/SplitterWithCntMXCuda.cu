// standard library includes
#include <iostream>
#include <math.h>

#include <gsl/gsl_cdf.h>

#include "SplitterWithCntMXCuda.cuh"
#include "kernels.cuh"

extern const float betaParam;
extern const float aParam;
extern const float bParam;
extern const double significanceThreshold;
extern const int minSeqLength;
extern const std::string symbols;

typedef unsigned char BYTE;

std::vector<std::string> SplitterWithCntMXCuda::split(int startIdx, int endIdx)
{
	double maxDivergence = -1.0;
	int seqLen = endIdx - startIdx + 1;
	int maxPos = -1;
	if (seqLen >= minSeqLength)
	{
		double * divs_A = computeDivergenceFromCountArray(countMatrixHolder.getCountArrayA(), startIdx, endIdx, seqLen);
		double * divs_C = computeDivergenceFromCountArray(countMatrixHolder.getCountArrayC(), startIdx, endIdx, seqLen);
		double * divs_G = computeDivergenceFromCountArray(countMatrixHolder.getCountArrayG(), startIdx, endIdx, seqLen);
		double * divs_T = computeDivergenceFromCountArray(countMatrixHolder.getCountArrayT(), startIdx, endIdx, seqLen);
		determineMaxDivAndPos(startIdx, seqLen, divs_A, divs_C, divs_G, divs_T, maxDivergence, maxPos);
	}
	double significance = computeSignificance(seqLen, maxDivergence);
	std::vector<std::string> result = checkSignificance(startIdx, endIdx, maxPos, significance);
	return result;
}

void SplitterWithCntMXCuda::determineMaxDivAndPos(int startIdx, int seqLen, double * divs_A, double * divs_C, double * divs_G, double * divs_T, double& maxDivergence, int& maxPos)
{
	double divergence;
	for (int pos = 1; pos < seqLen; pos++)
	{
		double divA = divs_A[pos];
		double divC = divs_C[pos];
		double divG = divs_G[pos];
		double divT = divs_T[pos];
		divergence = divA + divC + divG + divT;
		if (maxDivergence < divergence)
		{
			maxDivergence = divergence;
			maxPos = startIdx + pos;
		}
	}
}

double * SplitterWithCntMXCuda::computeDivergenceFromCountArray(BYTE * countArray, int startIdx, int endIdx, int seqLen)
{
	// CUDA grid launch parameters
	int nThreads = 256;
	dim3 nT(nThreads);
	dim3 nB((seqLen + nThreads - 1) / nThreads);
	if (nB.x > 65535)
	{
		std::stringstream errmsg;
		errmsg << "ERROR: Block is too large:\n";
		errmsg << nB.x << " blocks. Max is 65535.\n";
		throw std::runtime_error(errmsg.str());
	}

	int countArraySize = countMatrixHolder.getCountArraySize();
	int size = countArraySize * sizeof(BYTE);
	BYTE * d_countArray;
	cudaMalloc((void **)&d_countArray, size);
	cudaCheckError(cudaMemcpy(d_countArray, countArray, size, cudaMemcpyHostToDevice));

	double * d_divs;
	int sizeOfDivs = seqLen * sizeof(double);
	cudaMalloc((void **)&d_divs, sizeOfDivs);

	computePartDivergenceForPos << < nB, nT >> >(d_countArray, startIdx, endIdx, seqLen, d_divs);

	cudaCheckError(cudaPeekAtLastError());
	cudaCheckError(cudaDeviceSynchronize());

	double * divs = new double[seqLen];
	cudaCheckError(cudaMemcpy(divs, d_divs, sizeOfDivs, cudaMemcpyDeviceToHost));

	cudaFree(d_countArray);
	cudaFree(d_divs);

	return divs;
}

double SplitterWithCntMXCuda::computeSignificance(int N, double maxDivergence)
{
	double significance = 0.0;
	if (maxDivergence > 0.0)
	{
		double NEff = aParam * log(N) + bParam;
		significance = gsl_cdf_chisq_P(2 * N * log(2) * betaParam * maxDivergence, symbols.length() - 1);
		significance = pow(significance, NEff);
	}
	return significance;
}

std::vector<std::string> SplitterWithCntMXCuda::checkSignificance(int startIdx, int endIdx, int maxPos, double significance)
{
	std::vector<std::string> result;
	if ((significance > significanceThreshold) && (maxPos >= 0))
	{
		std::cout << "New splitpoint on pos " << maxPos << std::endl;
		std::vector<std::string>& resultPrefix = split(startIdx, maxPos - 1);
		std::vector<std::string>& resultPostfix = split(maxPos, endIdx);
		result.resize(resultPrefix.size() + resultPostfix.size());
		result = std::vector<std::string>(resultPrefix);
		result.insert(result.end(), resultPostfix.begin(), resultPostfix.end());
	}
	else
	{
		result.resize(1);
		std::string resultSeq = sequence.substr(startIdx, endIdx - startIdx + 1);
		result.push_back(resultSeq);
	}
	return result;
}
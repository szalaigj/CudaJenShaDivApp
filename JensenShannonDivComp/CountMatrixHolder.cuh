#ifndef COUNTMXHOLDER_HPP_
#define COUNTMXHOLDER_HPP_

#include <math.h>
#include <string>
#include <algorithm>

// std includes
#include <exception>
#include <sstream>

#include "kernels.cuh"

// local CudaJenShaDivApp includes
#include "CudaBasicIncludes.cuh"

typedef unsigned char BYTE;

extern const BYTE subSeqLength;

class ICountMatrixHolder
{
public:
	virtual int queryCount(char symbol, int i, int j) = 0;
	virtual double queryFrequency(char symbol, int i, int j) = 0;
};

class BasicCountMatrixHolder : public ICountMatrixHolder
{
public:
	const BYTE * getCountArrayA() const
	{
		return countArrayA;
	}

	const BYTE * getCountArrayC() const
	{
		return countArrayC;
	}

	const BYTE * getCountArrayG() const
	{
		return countArrayG;
	}

	const BYTE * getCountArrayT() const
	{
		return countArrayT;
	}

	virtual void innerSetupCountArrays(std::string& sequence, int seqLen, int subSeqNumber, int countArraySize) = 0;
	int transformIndex(int k, int i, int j);
	BYTE queryCountInSubSeq(char symbol, int k, int i, int j);
	int queryCount(char symbol, int i, int j);
	double queryFrequency(char symbol, int i, int j);

	void setupCountArrays(std::string& sequence)
	{
		int seqLen = sequence.length();
		int subSeqNumber = ceil((double)seqLen / (double)subSeqLength);
		int countArraySize = subSeqNumber * ((subSeqLength + 1) * subSeqLength / 2);
		countArrayA = new BYTE[countArraySize];
		countArrayC = new BYTE[countArraySize];
		countArrayG = new BYTE[countArraySize];
		countArrayT = new BYTE[countArraySize];
		innerSetupCountArrays(sequence, seqLen, subSeqNumber, countArraySize);
	}
protected:
	BYTE * countArrayA;
	BYTE * countArrayC;
	BYTE * countArrayG;
	BYTE * countArrayT;
};

class CountMatrixHolder : public BasicCountMatrixHolder
{
public:
	void innerSetupCountArrays(std::string& sequence, int seqLen, int subSeqNumber, int countArraySize);
	void count(std::string& currentInnerSubSeq, BYTE& a, BYTE& c, BYTE& g, BYTE& t);
};

class CountMatrixHolderCuda : public BasicCountMatrixHolder
{
public:
	void innerSetupCountArrays(std::string& sequence, int seqLen, int subSeqNumber, int countArraySize);
};
#endif /* COUNTMXHOLDER_HPP_ */
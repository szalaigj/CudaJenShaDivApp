#include "CountMatrixHolder.cuh"

extern const BYTE subSeqLength;

int BasicCountMatrixHolder::transformIndex(int k, int i, int j)
{
	return k * ((subSeqLength + 1) * subSeqLength / 2) + (j + 1) * j / 2 + i;
}

BYTE BasicCountMatrixHolder::queryCountInSubSeq(char symbol, int k, int i, int j)
{
	const BYTE * currentCountArray = nullptr;
	switch (symbol)
	{
	case 'A':
		currentCountArray = getCountArrayA();
		break;
	case 'C':
		currentCountArray = getCountArrayC();
		break;
	case 'G':
		currentCountArray = getCountArrayG();
		break;
	case 'T':
		currentCountArray = getCountArrayT();
		break;
	default:
		break;
	}
	BYTE result = currentCountArray[transformIndex(k, i, j)];
	return result;
}

int BasicCountMatrixHolder::queryCount(char symbol, int i, int j)
{
	int temp = i;
	if (i > j)
	{
		i = j;
		j = temp;
	}
	int I = i % subSeqLength;
	int J = j % subSeqLength;
	int Ki = (i - I) / subSeqLength;
	int Kj = (j - J) / subSeqLength;
	int result;
	if (Ki == Kj)
	{
		result = queryCountInSubSeq(symbol, Ki, I, J);
	}
	else if (Kj == Ki + 1)
	{
		result = queryCountInSubSeq(symbol, Ki, I, subSeqLength - 1);
		result += queryCountInSubSeq(symbol, Kj, 0, J);
	}
	else
	{
		result = queryCountInSubSeq(symbol, Ki, I, subSeqLength - 1);
		for (int k = Ki + 1; k <= Kj - 1; k++)
		{
			result += queryCountInSubSeq(symbol, k, 0, subSeqLength - 1);
		}
		result += queryCountInSubSeq(symbol, Kj, 0, J);
	}
	return result;
}

double BasicCountMatrixHolder::queryFrequency(char symbol, int i, int j)
{
	int temp = i;
	if (i > j)
	{
		i = j;
		j = temp;
	}
	double result;
	result = (double)queryCount(symbol, i, j) / (double)(j - i + 1);
	return result;
}

void CountMatrixHolder::innerSetupCountArrays(std::string& sequence, int seqLen, int subSeqNumber, int countArraySize)
{
	BYTE currentCount;
	std::string currentOuterSubSeq;
	std::string currentInnerSubSeq;
	for (int k = 0; k < subSeqNumber; k++)
	{
		int currentOuterStartIdx = k * subSeqLength;
		int currentSubSeqLength = std::min((int)subSeqLength, seqLen - currentOuterStartIdx + 1);
		currentOuterSubSeq = sequence.substr(currentOuterStartIdx, currentSubSeqLength);
		for (int j = 0; j < currentSubSeqLength; j++)
		{
			for (int i = 0; i <= j; i++)
			{
				currentInnerSubSeq = currentOuterSubSeq.substr(i, j - i + 1);
				BYTE a = 0;
				BYTE c = 0;
				BYTE g = 0;
				BYTE t = 0;
				count(currentInnerSubSeq, a, c, g, t);
				int currentIdxOfCountArrays = transformIndex(k, i, j);
				countArrayA[currentIdxOfCountArrays] = a;
				countArrayC[currentIdxOfCountArrays] = c;
				countArrayG[currentIdxOfCountArrays] = g;
				countArrayT[currentIdxOfCountArrays] = t;
				//printf("(k, j, i) : (%i, %i, %i); (a, c, g, t) : (%i, %i, %i, %i)\n", k, j, i, a, c, g, t);
			}
		}
	}
}

void CountMatrixHolder::count(std::string& currentInnerSubSeq, BYTE& a, BYTE& c, BYTE& g, BYTE& t)
{
	for (int idx = 0; idx < currentInnerSubSeq.length(); idx++)
	{
		char currentChr = currentInnerSubSeq[idx];
		switch (currentChr)
		{
		case 'A':
			a++;
			break;
		case 'C':
			c++;
			break;
		case 'G':
			g++;
			break;
		case 'T':
			t++;
			break;
		default:
			break;
		}
	}
}
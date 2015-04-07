#include "kernels.cuh"

extern __constant__ BYTE d_subSeqLength;
extern __constant__ int d_sequenceLength;
extern __constant__ int d_subSeqNumber;

__device__
void countKernel(char * tmp, int j, int i, BYTE& a, BYTE& c, BYTE& g, BYTE& t)
{
	for (int idx = i; idx <= j; idx++)
	{
		char currentChr = tmp[idx];
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

__device__
int transformIndex(int k, int i, int j)
{
	return k * ((d_subSeqLength + 1) * d_subSeqLength / 2) + (j + 1) * j / 2 + i;
}

__global__
void setupCountArraysKernel(char * sequence, BYTE * countArrayA, BYTE * countArrayC, BYTE * countArrayG, BYTE * countArrayT)
{
	__shared__ char tmp[255];
	int k = blockIdx.x;
	int j = threadIdx.x;
	int maximalCurrentSubSeqLength = min((int)d_subSeqLength, d_sequenceLength - k * d_subSeqLength + 1);
	
	if (j < maximalCurrentSubSeqLength)
	{
		tmp[j] = sequence[k * d_subSeqLength + j];
	}

	__syncthreads();

	if ((k < d_subSeqNumber) && (j < maximalCurrentSubSeqLength))
	{
		for (int i = 0; i <= j; i++)
		{
			BYTE a = 0;
			BYTE c = 0;
			BYTE g = 0;
			BYTE t = 0;
			
			countKernel(tmp, j, i, a, c, g, t);

			int currentIdxOfCountArrays = transformIndex(k, i, j);
			countArrayA[currentIdxOfCountArrays] = a;
			countArrayC[currentIdxOfCountArrays] = c;
			countArrayG[currentIdxOfCountArrays] = g;
			countArrayT[currentIdxOfCountArrays] = t;
			//printf("(k, j, i) : (%i, %i, %i); (a, c, g, t) : (%u, %u, %u, %u)\n", k, j, i, a, c, g, t);
		}

	}
}

__device__
int queryCount(BYTE * countArray, int i, int j)
{
	int temp = i;
	if (i > j)
	{
		i = j;
		j = temp;
	}
	int I = i % d_subSeqLength;
	int J = j % d_subSeqLength;
	int Ki = (i - I) / d_subSeqLength;
	int Kj = (j - J) / d_subSeqLength;
	int result;
	if (Ki == Kj)
	{
		result = countArray[transformIndex(Ki, I, J)]; //queryCountInSubSeq(symbol, Ki, I, J);
	}
	else if (Kj == Ki + 1)
	{
		result = countArray[transformIndex(Ki, I, d_subSeqLength - 1)];
		result += countArray[transformIndex(Kj, 0, J)];
	}
	else
	{
		result = countArray[transformIndex(Ki, I, d_subSeqLength - 1)];
		for (int k = Ki + 1; k <= Kj - 1; k++)
		{
			result += countArray[transformIndex(k, 0, d_subSeqLength - 1)];
		}
		result += countArray[transformIndex(Kj, 0, J)];
	}
	return result;
}

__device__
double queryFrequency(BYTE * countArray, int i, int j)
{
	int temp = i;
	if (i > j)
	{
		i = j;
		j = temp;
	}
	double result;
	result = (double)queryCount(countArray, i, j) / (double)(j - i + 1);
	return result;
}

__device__
double computeEntropy(double frequency)
{
	double entropy = 0.0;
	if (frequency != 0.0)
	{
		entropy -= frequency * log(frequency);
	}
	return entropy;
}

__device__
double computeDivergence(double frequencyPrefix, double frequencyPostfix, double weightPrefix, double weightPostfix)
{
	double weightedFrequency1 = frequencyPrefix * weightPrefix;
	double weightedFrequency2 = frequencyPostfix * weightPostfix;
	double sumOfFrequencies = weightedFrequency1 + weightedFrequency2;
	double result = computeEntropy(sumOfFrequencies)
		- weightPrefix * computeEntropy(frequencyPrefix)
		- weightPostfix * computeEntropy(frequencyPostfix);
	return result;
}

__global__
void computePartDivergenceForPos(BYTE * countArray, int startIdx, int endIdx, int seqLen, double * divs)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	if ((index > 0) && (index < seqLen))
	{
		double frequencyPrefix = queryFrequency(countArray, startIdx, startIdx + index - 1);
		double frequencyPostfix = queryFrequency(countArray, startIdx + index, endIdx);
		double weightPrefix = (double)(index) / (double)(seqLen);
		double weightPostfix = (double)(seqLen - index) / (double)(seqLen);
		divs[index] = computeDivergence(frequencyPrefix, frequencyPostfix, weightPrefix, weightPostfix);
	}
}
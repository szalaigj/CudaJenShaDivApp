#include "kernels.cuh"

extern __constant__ BYTE d_subSeqLength;
extern __constant__ int d_sequenceLength;
extern __constant__ int d_subSeqNumber;
//extern __device__ char * d_sequence;

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
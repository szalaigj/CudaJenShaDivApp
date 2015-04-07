#ifndef _KERNELS_CUH__
#define _KERNELS_CUH__

// local CudaJenShaDivApp includes
#include "CudaBasicIncludes.cuh"

typedef unsigned char BYTE;

extern __constant__ BYTE d_subSeqLength;
extern __constant__ int d_sequenceLength;
extern __constant__ int d_subSeqNumber;
extern __device__ char * d_sequence;

__device__
void countKernel(char * tmp, int j, int i, BYTE& a, BYTE& c, BYTE& g, BYTE& t);

__device__
int transformIndex(int k, int i, int j);

__global__
void setupCountArraysKernel(char * sequence, BYTE * countArrayA, BYTE * countArrayC, BYTE * countArrayG, BYTE * countArrayT);

__device__
int queryCount(BYTE * countArray, int i, int j);

__device__
double queryFrequency(BYTE * countArray, int i, int j);

__device__
double computeEntropy(double frequency);

__device__
double computeDivergence(double frequencyPrefix, double frequencyPostfix, double weightPrefix, double weightPostfix);

__global__
void computePartDivergenceForPos(BYTE * countArray, int startIdx, int endIdx, int seqLen, double * divs);

#endif /* _KERNELS_CUH__ */
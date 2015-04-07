#include "CountMatrixHolder.cuh"

typedef unsigned char BYTE;

extern const BYTE subSeqLength;

__constant__ BYTE d_subSeqLength;
__constant__ int d_sequenceLength;
__constant__ int d_subSeqNumber;

void CountMatrixHolderCuda::innerSetupCountArrays(std::string& sequence, int seqLen, int subSeqNumber, int countArraySize)
{
	// CUDA grid launch parameters
	int nThreads = 256;
	dim3 nT(nThreads);
	dim3 nB((seqLen + subSeqLength - 1) / subSeqLength);
	if (nB.x > 65535)
	{
		std::stringstream errmsg;
		errmsg << "ERROR: Block is too large:\n";
		errmsg << nB.x << " blocks. Max is 65535.\n";
		throw std::runtime_error(errmsg.str());
	}

	cudaCheckError(cudaMemcpyToSymbol(d_subSeqLength, &subSeqLength, sizeof(BYTE)));
	cudaCheckError(cudaMemcpyToSymbol(d_sequenceLength, &seqLen, sizeof(int)));
	cudaCheckError(cudaMemcpyToSymbol(d_subSeqNumber, &subSeqNumber, sizeof(int)));
	//cudaCheckError(cudaMemcpyToSymbol(d_sequence, sequence.c_str(), (seqLen + 1) * sizeof(*(sequence.c_str()))));
	char * d_sequence;
	int sizeOfSeq = (seqLen + 1) * sizeof(char);
	cudaMalloc((void **)&d_sequence, sizeOfSeq);
	char * seqChars = (char *)sequence.c_str();
	cudaCheckError(cudaMemcpy(d_sequence, seqChars, sizeOfSeq, cudaMemcpyHostToDevice));

	BYTE * d_countArrayA;
	BYTE * d_countArrayC;
	BYTE * d_countArrayG;
	BYTE * d_countArrayT;

	int size = countArraySize * sizeof(BYTE);
	cudaMalloc((void **)&d_countArrayA, size);
	cudaMalloc((void **)&d_countArrayC, size);
	cudaMalloc((void **)&d_countArrayG, size);
	cudaMalloc((void **)&d_countArrayT, size);

	setupCountArraysKernel << < nB, nT >> >(d_sequence, d_countArrayA, d_countArrayC, d_countArrayG, d_countArrayT);

	cudaCheckError(cudaPeekAtLastError());
	cudaCheckError(cudaDeviceSynchronize());

	cudaCheckError(cudaMemcpy(countArrayA, d_countArrayA, size, cudaMemcpyDeviceToHost));
	cudaCheckError(cudaMemcpy(countArrayC, d_countArrayC, size, cudaMemcpyDeviceToHost));
	cudaCheckError(cudaMemcpy(countArrayG, d_countArrayG, size, cudaMemcpyDeviceToHost));
	cudaCheckError(cudaMemcpy(countArrayT, d_countArrayT, size, cudaMemcpyDeviceToHost));

	cudaFree(d_sequence);
	cudaFree(d_countArrayA);
	cudaFree(d_countArrayC);
	cudaFree(d_countArrayG);
	cudaFree(d_countArrayT);
}
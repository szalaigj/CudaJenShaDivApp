#ifndef _CUDA_BASIC_INCLUDES__
#define _CUDA_BASIC_INCLUDES__

// Cuda Includes
#include "cuda_runtime.h"
#include "cuda_runtime_api.h"
#include "device_launch_parameters.h"
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>

// standard includes
#include <fstream>

inline
void cudaAssert(cudaError_t code, char * file, int line)
{
	if (code != 0) {
		fprintf(stderr, "An error '%s' occured in file '%s' at line %d\n",
			cudaGetErrorString(code), file, line);
		exit(1);
	}
}
#define cudaCheckError(ans) { cudaAssert((ans), __FILE__, __LINE__); }

#endif // _CUDA_BASIC_INCLUDES__ //
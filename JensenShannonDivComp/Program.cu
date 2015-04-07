/*
* This implementation is based on the article:
*	Grosse, I., Bernaola-Galvan, P., Carpena, P., Roman-Roldan, R., Oliver, J., & Stanley, H. E. (2002).
*	Analysis of symbolic sequences using the Jensen-Shannon divergence.
*	Physical Review E, 65(4), 041905
*/

// standard library includes
#include <iostream>
#include <time.h>

// local CudaJenShaDivApp includes
#include "CudaBasicIncludes.cuh"

#include "../Utils/ShannonEntropyComputer.cpp"
#include "../Utils/FrequencyComputer.cpp"
#include "../Utils/JenShaDivComputer.cpp"
#include "../DataHandlingUtil/LoadData.cpp"
#include "../DataHandlingUtil/SubSeqWriter.cpp"
#include "../Utils/Splitter.cpp"
#include "../Utils/ChiSquaredCDFComputer.cpp"
#include "../Utils/UIncGamma.cpp"

#include "SplitterWithCntMX.hpp"
#include "CountMatrixHolder.cuh"
#include "SplitterWithCntMXCuda.cuh"

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

typedef unsigned char BYTE;
const BYTE subSeqLength = 255;

__constant__ double usedParams[3];

int main(int argc, char** argv)
{
	clock_t begin, end;
	double elapsed_secs;
	std::string seq;
	//seq = "GAATTCATTCACCATTATTCTTTTATAATATTGCTATTTTATTATTCTTGATCATTGACATTACCATATAAATGGTAGAATCAGCTTGAAAATTTCTACCAAAATGCCGCTTGGAATTTTTATTAGAATTGCATTGGATCTGGAGATCAATTTACGAAGAACTGACTCTTTAAACATAACAACTCTTCTGATCCATGACAAGGTTTATCTCCCCACTAATTTAGTTCTTTCATAATTTCTCAAAGCAATTTTTTGTAGTTTTTGGTGTACTGGCCTTACATAAATTTTGTTGACTTTCCTTTTTTTTTTTTTTTTTTTTTTTTTTTTGAGACAGAGTCTTGCTCTGTTACCCAGGCTGGAGTGCAGTGATGCGATCTCGGCTCACTGTAACCTCTGCCTCCCAGGTTCAAGTGATCCTCCTGCCTCAGCCTCCCAAGTAGCTGGACTACAGGCACATGCCACCACGCCTAGCTAATTTTTTGTATTTTTAGTAGAGATGGGGTTTCACTGTGTTAGGCGTGATGGTCTCCATCTCCTGACCTCATGATCTGCCCACCTCGGCTTCCCAAAGTGCTGAGATTACAAGTGTGAGCCATGGTGCCCAGTCTGTTGAATTTATTTATAAGCACAACATGTATTTAGATGTTACTTTAAATGAAATTGTATTCTTATTTCATTTTCCAAATGCTCATTGCTAATATACAGAAATACAAAAGACCACTTACATTGAGAGCTTACATTCTGCAACACTACCAAACTCACTGATTAGTTCTGGTAGATTTTTGTAGATTTCTAGCATTGTTAACAAACACAGTCATTATCTGTGAATAAAGACAGCTTCAATTCTTTCAATACTTTATTAATTTTTCTTACTTTATTGCATTGACTTAGAGCTCTAGTATAATGCTGAATTAAAAGAGTAACAACAGGTATTCTACTTTTTTCTCTGATTTAATAGAAAAGCATTCAATCTTATGCCATTTAATATAATGTTACCTGTGGGTTCTTCAAATCTGCCCTTAATGGGGTTGGAAGTGTTGCCTTCTGTTCTCATCATGCTGAGCATTTTCTGGGGTTTGTTTTTATAAATCATGAAAAAAGTTTTCAATTTTGCCAAATGCTTTTACTGTGTATGACAAGGTAATCATACGGTTTTTCTCTTTTGCCCTGATAATATATAAAATGATATTTTTTAAATATAAAAAAGAGGCCGGGCATGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCTGAGGCGGGCGGATCACACTTGTGGCAGCATTGAAGGCTTCACTCTTCCCCAAGGG";
	BaseDataLoader loader;
	ShannonEntropyComputer entropyComputer;
	FrequencyComputer frequencyComputer(symbols);
	JenShaDivComputer jenShaDivComputer(entropyComputer);
	UIncGamma uIncGamma;
	ChiSquaredCDFComputer chiSquaredCDFComputer(uIncGamma);
	Splitter splitter(frequencyComputer, jenShaDivComputer, chiSquaredCDFComputer);
	BaseSubSeqWriter subSeqWriter;
	
	std::string filename(argv[1]);
	seq = loader.loadData(filename);

	std::string outputPath(argv[2]);

	std::vector<std::string> outputPerformance;
	
	std::string outputPerfStr;
	std::ostringstream s;

	begin = clock();
	CountMatrixHolderCuda countMatrixHolderCuda0;
	countMatrixHolderCuda0.setupCountArrays(seq);
	end = clock();
	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	std::cout << "Part elapsed time (in min) of count matrix - splitter CUDA method" << " " << (elapsed_secs / 60.0) << std::endl;
	s << "Part elapsed time (in min) of count matrix - splitter CUDA method" << " " << (elapsed_secs / 60.0) << std::endl;
	outputPerfStr = s.str();
	outputPerformance.push_back(outputPerfStr);
	s.str("");

	SplitterWithCntMXCuda splitterWithCntMXCuda(seq, countMatrixHolderCuda0);
	std::string seqOutputFilename0(outputPath + "seq_output0.txt");
	std::vector<std::string> subsequences0 = splitterWithCntMXCuda.split(0, seq.length() - 1);
	subSeqWriter.writeData(seqOutputFilename0, &subsequences0);
	end = clock();
	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	std::cout << "Elapsed time (in min) of count matrix - splitter CUDA method" << " " << (elapsed_secs / 60.0) << std::endl;
	s << "Elapsed time (in min) of count matrix - splitter CUDA method" << " " << (elapsed_secs / 60.0) << std::endl;
	outputPerfStr = s.str();
	outputPerformance.push_back(outputPerfStr);
	s.str("");


	begin = clock();
	std::vector<std::string> subsequences = splitter.split(seq);
	std::string seqOutputFilename1(outputPath + "seq_output1.txt");
	subSeqWriter.writeData(seqOutputFilename1, &subsequences);
	end = clock();
	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	std::cout << "Elapsed time (in min) of original method" << " " << (elapsed_secs / 60.0) << std::endl;
	s << "Elapsed time (in min) of original method" << " " << (elapsed_secs / 60.0) << std::endl;
	outputPerfStr = s.str();
	outputPerformance.push_back(outputPerfStr);
	s.str("");

	begin = clock();
	CountMatrixHolder countMatrixHolder;
	countMatrixHolder.setupCountArrays(seq);

	//int aCount1 = countMatrixHolder.queryCount('A', 0, 547495);//163078
	//int aCount1 = countMatrixHolder.queryCount('A', 0, 1296);//348
	//int cCount1 = countMatrixHolder.queryCount('C', 0, 1310);//252

	end = clock();
	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	std::cout << "Part elapsed time (in min) of count matrix method" << " " << (elapsed_secs / 60.0) << std::endl;
	s << "Part elapsed time (in min) of count matrix method" << " " << (elapsed_secs / 60.0) << std::endl;
	outputPerfStr = s.str();
	outputPerformance.push_back(outputPerfStr);
	s.str("");

	SplitterWithCntMX splitterWithCntMX(seq, countMatrixHolder, jenShaDivComputer);
	std::string seqOutputFilename2(outputPath + "seq_output2.txt");
	std::vector<std::string> subsequences2 = splitterWithCntMX.split(0, seq.length() - 1);
	subSeqWriter.writeData(seqOutputFilename2, &subsequences2);
	end = clock();
	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	std::cout << "Elapsed time (in min) of count matrix method" << " " << (elapsed_secs / 60.0) << std::endl;
	s << "Elapsed time (in min) of count matrix method" << " " << (elapsed_secs / 60.0) << std::endl;
	outputPerfStr = s.str();
	outputPerformance.push_back(outputPerfStr);
	s.str("");

	begin = clock();
	CountMatrixHolderCuda countMatrixHolderCuda;
	countMatrixHolderCuda.setupCountArrays(seq);

	end = clock();
	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	std::cout << "Part elapsed time (in min) of count matrix CUDA method" << " " << (elapsed_secs / 60.0) << std::endl;
	s << "Part elapsed time (in min) of count matrix CUDA method" << " " << (elapsed_secs / 60.0) << std::endl;
	outputPerfStr = s.str();
	outputPerformance.push_back(outputPerfStr);
	s.str("");

	SplitterWithCntMX splitterWithCntMX2(seq, countMatrixHolderCuda, jenShaDivComputer);
	std::string seqOutputFilename3(outputPath + "seq_output3.txt");
	std::vector<std::string> subsequences3 = splitterWithCntMX2.split(0, seq.length() - 1);
	subSeqWriter.writeData(seqOutputFilename3, &subsequences3);
	end = clock();
	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	std::cout << "Elapsed time (in min) of count matrix CUDA method" << " " << (elapsed_secs / 60.0) << std::endl;
	s << "Elapsed time (in min) of count matrix CUDA method" << " " << (elapsed_secs / 60.0) << std::endl;
	outputPerfStr = s.str();
	outputPerformance.push_back(outputPerfStr);
	s.str("");

	std::string perfOutputFilename(outputPath + "performance.txt");
	subSeqWriter.writeData(perfOutputFilename, &outputPerformance);

	std::cout << "Press any character and press enter to continue..." << std::endl;
	char chr;
	std::cin >> chr;
}
#ifndef SPLITTER_WITH_CNT_MX_HPP_
#define SPLITTER_WITH_CNT_MX_HPP_

#include <vector>

#include "../Utils/JenShaDivComputer.hpp"
#include "CountMatrixHolder.cuh"

extern const std::string symbols;

class ISplitterWithCntMX
{
public:
	virtual std::vector<std::string> split(int startIdx, int endIdx) = 0;
};


class SplitterWithCntMX : public ISplitterWithCntMX
{
public:
	SplitterWithCntMX(std::string& sequence, ICountMatrixHolder& countMatrixHolder, JenShaDivComputer& jenShaDivComputer) :
		sequence(sequence), countMatrixHolder(countMatrixHolder), jenShaDivComputer(jenShaDivComputer)
	{
	}

	ICountMatrixHolder& getCountMatrixHolder()
	{
		return countMatrixHolder;
	}

	JenShaDivComputer& getJenShaDivComputer()
	{
		return jenShaDivComputer;
	}

	std::vector<std::string> split(int startIdx, int endIdx);
	void computeDivergenceForPos(int, int, int, double&, int&);
	double computeSignificance(int, double);
	std::vector<std::string> checkSignificance(int, int, int, double);

private:
	std::string& sequence;
	ICountMatrixHolder& countMatrixHolder;
	JenShaDivComputer& jenShaDivComputer;
};

#endif /* SPLITTER_WITH_CNT_MX_HPP_ */
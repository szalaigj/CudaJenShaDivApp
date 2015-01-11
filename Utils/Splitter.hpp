#ifndef SPLITTER_HPP_
#define SPLITTER_HPP_

#include <map>
#include <vector>

#include "FrequencyComputer.hpp"
#include "JenShaDivComputer.hpp"
#include "ChiSquaredCDFComputer.hpp"

extern const std::string symbols;

class ISplitter
{
	public:
		virtual std::vector<std::string> split(std::string& sequence) = 0;
		virtual void computeDivergenceForPos(std::string& sequence, double& maxDivergence,
			std::string& seqPrefixForMax, std::string& seqPostfixForMax, int pos,
			std::map <char, long>& chrCountsPrefix, std::map <char, long>& chrCountsPostfix) = 0;
		virtual double computeSignificance(int N, double maxDivergence) = 0;
		virtual std::vector<std::string> checkSignificance(std::string& sequence, std::string& seqPrefixForMax,
			std::string& seqPostfixForMax, double significance) = 0;
};

class Splitter : public ISplitter
{
	public:
		Splitter(FrequencyComputer& frequencyComputer, JenShaDivComputer& jenShaDivComputer,
			ChiSquaredCDFComputer& chiSquaredCDFComputer) :
			frequencyComputer(frequencyComputer), jenShaDivComputer(jenShaDivComputer),
			chiSquaredCDFComputer(chiSquaredCDFComputer)
		{
		}

		FrequencyComputer& getFrequencyComputer()
		{
			return frequencyComputer;
		}

		JenShaDivComputer& getJenShaDivComputer()
		{
			return jenShaDivComputer;
		}

		ChiSquaredCDFComputer& getChiSquaredCDFComputer()
		{
			return chiSquaredCDFComputer;
		}

		std::vector<std::string> split(std::string&);
		void computeDivergenceForPos(std::string&, double&, std::string&, std::string&, int,
			std::map <char, long>&, std::map <char, long>&);
		double computeSignificance(int, double);
		std::vector<std::string> checkSignificance(std::string&, std::string&, std::string&, double);

	private:
		FrequencyComputer& frequencyComputer;
		JenShaDivComputer& jenShaDivComputer;
		ChiSquaredCDFComputer& chiSquaredCDFComputer;
};

#endif /* SPLITTER_HPP_ */
#include "FrequencyComputer.hpp"

double * FrequencyComputer::computeFrequency(std::string& sequence)
{
	std::string symbols = getSymbols();
	double * frequencies = new double[symbols.length()];
	std::map <char, long> chrCounts;
	initChrCounts(chrCounts);
	countCharsInSeq(sequence, chrCounts);
	determineFrequencies(frequencies, chrCounts, sequence.length());
	return frequencies;
}

void FrequencyComputer::initChrCounts(std::map <char, long>& chrCounts)
{
	std::string symbols = getSymbols();
	for (int idx = 0; idx < symbols.length(); idx++)
	{
		chrCounts[symbols[idx]] = 0L;
	}
}

void FrequencyComputer::countCharsInSeq(std::string& sequence, std::map <char, long>& chrCounts)
{
	for (int idx = 0; idx < sequence.length(); idx++)
	{
		char currentChr = sequence[idx];
		if (chrCounts.find(currentChr) == chrCounts.end())
		{
			throwErrorMsg(currentChr);
		}
		else
		{
			chrCounts[currentChr]++;
		}
	}
}

void FrequencyComputer::increaseCountChar(char chrInSeq, std::map <char, long>& chrCounts)
{
	if (chrCounts.find(chrInSeq) == chrCounts.end())
	{
		throwErrorMsg(chrInSeq);
	}
	else
	{
		chrCounts[chrInSeq]++;
	}
}

void FrequencyComputer::decreaseCountChar(char chrInSeq, std::map <char, long>& chrCounts)
{
	if (chrCounts.find(chrInSeq) == chrCounts.end())
	{
		throwErrorMsg(chrInSeq);
	}
	else
	{
		chrCounts[chrInSeq]--;
	}
}

double * FrequencyComputer::determineFrequencies(std::map <char, long>& chrCounts, int sequenceLength)
{
	std::string symbols = getSymbols();
	double * frequencies = new double[symbols.length()];
	for (int idx = 0; idx < symbols.length(); idx++)
	{
		frequencies[idx] = (double)chrCounts[symbols[idx]] / (double)sequenceLength;
	}
	return frequencies;
}